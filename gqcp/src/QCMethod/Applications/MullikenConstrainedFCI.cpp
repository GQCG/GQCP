// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2019  the GQCG developers
//
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
//
#include "QCMethod/Applications/MullikenConstrainedFCI.hpp"

#include "Basis/transform.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Mathematical/Optimization/Eigenproblem/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Processing/Properties/expectation_values.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"

#include <algorithm>
#include <chrono>


namespace GQCP {
namespace QCMethod {

/*
 *  PRIVATE METHODS
 */

/**
 *  Store the solutions from a solve
 *  
 *  @param eigenpairs           the eigenpairs from the CI solver
 *  @param multiplier           the Lagrangian multiplier associated with the solution
 *  @param sz_multiplier        a given multiplier for the atomic Sz constraint
 */ 
void MullikenConstrainedFCI::parseSolution(const std::vector<Eigenpair>& eigenpairs, const double multiplier, const double sz_multiplier) {

    // Initialize the result vectors to zero
    if (this->energy.size() != eigenpairs.size()) {
        this->energy = std::vector<double>(eigenpairs.size());
        this->population = std::vector<double>(eigenpairs.size());
        this->lambda = std::vector<double>(eigenpairs.size());
        this->lambda_sz = std::vector<double>(eigenpairs.size());
        this->entropy = std::vector<double>(eigenpairs.size());
        this->sz = std::vector<double>(eigenpairs.size());

        if (molecule.numberOfAtoms() == 2) {
            this->A_fragment_energy = std::vector<double>(eigenpairs.size());
            this->A_fragment_self_energy = std::vector<double>(eigenpairs.size());
            this->B_fragment_energy = std::vector<double>(eigenpairs.size());
            this->B_fragment_self_energy = std::vector<double>(eigenpairs.size());
            this->interaction_energy = std::vector<double>(eigenpairs.size());
        }

        this->eigenvector = std::vector<VectorX<double>>(eigenpairs.size());
    }

    // Fill in the results
    double internuclear_repulsion_energy = Operator::NuclearRepulsion(this->molecule).value();

    for (size_t i = 0; i < eigenpairs.size(); i++) {

        const auto& pair = eigenpairs[i];
        const auto& fci_coefficients = pair.get_eigenvector();
        double fci_energy = pair.get_eigenvalue();
        this->rdm_calculator.set_coefficients(fci_coefficients);
        const auto rdms = this->rdm_calculator.calculate1RDMs();
        OneRDM<double> D = rdms.one_rdm;
        OneRDM<double> D_s = rdms.one_rdm_aa - rdms.one_rdm_bb;
        TwoRDM<double> d = this->rdm_calculator.calculate2RDMs().two_rdm;

        double population = mulliken_operator.calculateExpectationValue(D)(0);
        double sz = sq_sz_operator.calculateExpectationValue(D_s)(0);
        LinearExpansion linear_expansion (fock_space, fci_coefficients);

        this->energy[i] = pair.get_eigenvalue() + internuclear_repulsion_energy + multiplier * population + sz_multiplier * sz;
        this->population[i] = population;
        this->sz[i] = sz;
        this->lambda[i] = multiplier;
        this->lambda_sz[i] = sz_multiplier;
        this->entropy[i] = linear_expansion.calculateShannonEntropy();

        if (molecule.numberOfAtoms() == 2) {
            // Transform the RDMs to the atomic orbital basis
            D.basisTransformInPlace(this->spinor_basis.coefficientMatrix().adjoint());
            d.basisTransformInPlace(this->spinor_basis.coefficientMatrix().adjoint());

            this->A_fragment_energy[i] = adp.get_atomic_parameters()[0].calculateExpectationValue(D, d) + internuclear_repulsion_energy / 2;
            this->A_fragment_self_energy[i] = adp.get_net_atomic_parameters()[0].calculateExpectationValue(D, d);
            this->B_fragment_energy[i] = adp.get_atomic_parameters()[1].calculateExpectationValue(D, d) + internuclear_repulsion_energy / 2;
            this->B_fragment_self_energy[i] = adp.get_net_atomic_parameters()[1].calculateExpectationValue(D, d);
            this->interaction_energy[i] = adp.get_interaction_parameters()[0].calculateExpectationValue(D, d) + internuclear_repulsion_energy;
        }

        this->eigenvector[i] = fci_coefficients;
    }
}


/**
 *  Throws an error if no solution is available
 *  
 *  @param function_name            name of the function that should throw the error
 */
void MullikenConstrainedFCI::checkAvailableSolutions(const std::string& function_name) const {
    if (!are_solutions_available) {
        throw std::runtime_error("MullikenConstrainedFCI::" + function_name + "(): The method hasn't been solved yet");
    }
}

/**
 *  Throws an error if the molecule is not diatomic
 *  
 *  @param function_name            name of the function that should throw the error
 */
void MullikenConstrainedFCI::checkDiatomicMolecule(const std::string& function_name) const {
    if (molecule.numberOfAtoms() != 2) {
        throw std::runtime_error("MullikenConstrainedFCI::" + function_name + "(): This property is only available for diatomic molecules");
    }
}

/*
 * CONSTRUCTORS
 */

/**
 *  @param molecule                 the molecule that will be solved for
 *  @param basis_set                the basisset that should be used
 *  @param basis_targets            the targeted basis functions for the constraint
 *  @param frozencores              the amount of frozen cores for the FCI calculation
 */
MullikenConstrainedFCI::MullikenConstrainedFCI(const Molecule& molecule, const std::string& basis_set, const std::vector<size_t>& basis_targets, const size_t frozencores) : 
        basis_targets (basis_targets),
        molecule (molecule),
        spinor_basis (RSpinorBasis<double, GTOShell>(molecule, basis_set)),
        uspinor_basis (USpinorBasis<double, GTOShell>(molecule, basis_set)),
        sq_hamiltonian (SQHamiltonian<double>::Molecular(this->spinor_basis, molecule)),  // in AO basis
        usq_hamiltonian (USQHamiltonian<double>::Molecular(this->uspinor_basis, molecule)),  // in AO basis
        basis_set (basis_set)
{

    if ((molecule.numberOfElectrons() % 2) > 0) {
        throw std::runtime_error("MullikenConstrainedFCI::MullikenConstrainedFCI(): This module is not available for an odd number of electrons");
    }

    auto K = this->spinor_basis.simpleDimension();
    auto N_P = this->molecule.numberOfElectrons()/2;

    try {
        // Try the foward approach of solving the RHF equations
        auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), this->sq_hamiltonian, this->spinor_basis.overlap().parameters());
        auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS(6, 6, 1.0e-12, 500);
        const GQCP::DiagonalRHFFockMatrixObjective<double> objective (this->sq_hamiltonian);
        const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment).groundStateParameters();

        basisTransform(this->spinor_basis, this->sq_hamiltonian, rhf_parameters.coefficientMatrix());

    } catch (const std::exception& e) {

        // If the DIIS does not converge, attempt to solve the RHF for the individuals atoms (if diatomic) and recombine the solutions to create a total canonical matrix
        // Starting from this new basis we re-attempt the regular DIIS
        // If all else fails perform Lowdin orthonormalization
        if (molecule.numberOfAtoms() == 2) {
            try {
                const std::vector<Nucleus>& atoms = molecule.nuclearFramework().nucleiAsVector();
                int charge = - molecule.numberOfElectrons() + molecule.nuclearFramework().totalNucleicCharge();
                Molecule mol_fraction1(std::vector<Nucleus>{atoms[0]}, charge);
                Molecule mol_fraction2(std::vector<Nucleus>{atoms[1]}, 0);

                RSpinorBasis<double, GTOShell> spinor_basis1 (mol_fraction1, basis_set);
                RSpinorBasis<double, GTOShell> spinor_basis2 (mol_fraction2, basis_set);

                auto ham_par1 = SQHamiltonian<double>::Molecular(spinor_basis1, mol_fraction1);  // in AO basis
                auto ham_par2 = SQHamiltonian<double>::Molecular(spinor_basis2, mol_fraction2);  // in AO basis

                // Perform DIIS RHF for individual fractions
                auto rhf_environment1 = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(mol_fraction1.numberOfElectrons(), ham_par1, spinor_basis1.overlap().parameters());
                auto diis_rhf_scf_solver1 = GQCP::RHFSCFSolver<double>::DIIS(6, 6, 1.0e-12, 500);
                const GQCP::DiagonalRHFFockMatrixObjective<double> objective1 (ham_par1);
                const auto rhf_parameters1 = GQCP::QCMethod::RHF<double>().optimize(objective1, diis_rhf_scf_solver1, rhf_environment1).groundStateParameters();

                auto rhf_environment2 = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(mol_fraction2.numberOfElectrons(), ham_par2, spinor_basis2.overlap().parameters());
                auto diis_rhf_scf_solver2 = GQCP::RHFSCFSolver<double>::DIIS(6, 6, 1.0e-12, 500);
                const GQCP::DiagonalRHFFockMatrixObjective<double> objective2 (ham_par2);
                const auto rhf_parameters2 = GQCP::QCMethod::RHF<double>().optimize(objective2, diis_rhf_scf_solver2, rhf_environment2).groundStateParameters();

                // Retrieve transformation from the solutions and transform the Hamiltonian
                size_t K1 = ham_par1.dimension();
                size_t K2 = ham_par2.dimension();

                // Recombine canonical matrices
                TransformationMatrix<double> T = Eigen::MatrixXd::Zero(K, K);
                T.topLeftCorner(K1, K1) += rhf_parameters1.coefficientMatrix();
                T.bottomRightCorner(K2, K2) += rhf_parameters2.coefficientMatrix();
                basisTransform(this->spinor_basis, this->sq_hamiltonian, T);


                // Attempt the DIIS for this basis
                try {
                    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), this->sq_hamiltonian, this->spinor_basis.overlap().parameters());
                    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS(6, 6, 1.0e-12, 500);
                    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (this->sq_hamiltonian);
                    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment).groundStateParameters();

                    basisTransform(this->spinor_basis, this->sq_hamiltonian, rhf_parameters.coefficientMatrix());


                } catch (const std::exception& e) {
                    const auto T = this->spinor_basis.lowdinOrthonormalizationMatrix();
                    basisTransform(this->spinor_basis, this->sq_hamiltonian, T);
                }


            } catch (const std::exception& e) {
                const auto T = this->spinor_basis.lowdinOrthonormalizationMatrix();
                basisTransform(this->spinor_basis, this->sq_hamiltonian, T);
            }

        } else {
            const auto T = this->spinor_basis.lowdinOrthonormalizationMatrix();
            basisTransform(this->spinor_basis, this->sq_hamiltonian, T);
        }
    }

    basisTransform(this->uspinor_basis, this->usq_hamiltonian, this->spinor_basis.coefficientMatrix());
    this->fock_space = SpinResolvedFrozenONVBasis(K, N_P, N_P, frozencores);
    this->mulliken_operator = this->spinor_basis.calculateMullikenOperator(basis_targets);
    this->sq_sz_operator = this->uspinor_basis.calculateAtomicSpinZ(basis_targets, SpinComponent::ALPHA);


    // Atomic Decomposition is only available for diatomic molecules
    if (molecule.numberOfAtoms() == 2) {
        this->adp = AtomicDecompositionParameters::Nuclear(molecule, basis_set);
    }

    this->rdm_calculator = RDMCalculator(fock_space);
}


/*
 * PUBLIC METHODS
 */

/**
 *  Solve the eigenvalue problem for a multiplier with the davidson algorithm
 *  
 *  @param multiplier           a given multiplier
 *  @param guess                supply a davidson guess
 *  @param sz_multiplier        a given multiplier for the atomic Sz constraint
 */
void MullikenConstrainedFCI::solveMullikenDavidson(const double multiplier, const VectorX<double>& guess, const double sz_multiplier) {

    auto start_time = std::chrono::high_resolution_clock::now();

    auto constrained_ham_par = this->usq_hamiltonian.constrain(this->mulliken_operator, multiplier, SpinComponent::ALPHA);
    constrained_ham_par = constrained_ham_par.constrain(this->mulliken_operator, multiplier, SpinComponent::BETA);
    constrained_ham_par = constrained_ham_par.constrain(this->sq_sz_operator, sz_multiplier, SpinComponent::ALPHA);
    constrained_ham_par = constrained_ham_par.constrain(this->sq_sz_operator, -sz_multiplier, SpinComponent::BETA);

    // Davidson solver
    DavidsonSolverOptions solver_options(guess);
    solver_options.convergence_threshold = this->convergence_threshold;
    solver_options.correction_threshold = this->correction_threshold;
    solver_options.maximum_subspace_dimension = this->maximum_subspace_dimension;
    solver_options.collapsed_subspace_dimension = this->collapsed_subspace_dimension;
    solver_options.maximum_number_of_iterations = this->maximum_number_of_iterations;

    VectorX<double> dia = this->fock_space.evaluateOperatorDiagonal(constrained_ham_par);
    VectorFunction<double> matrixVectorProduct = [this, &constrained_ham_par, &dia](const GQCP::VectorX<double>& x) { return this->fock_space.evaluateOperatorMatrixVectorProduct(constrained_ham_par, x, dia); };
    DavidsonSolver solver (matrixVectorProduct, dia, solver_options);

    try {
        solver.solve();
    } catch (const std::exception& e) {
        std::cout << e.what() << "multiplier: " << multiplier;
        return;
    }

    this->parseSolution(solver.get_eigenpairs(), multiplier, sz_multiplier);
    this->are_solutions_available = true;

    auto stop_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = stop_time - start_time;  // in nanoseconds
    this->solve_time = static_cast<double>(elapsed_time.count() / 1e9);  // in seconds
}

/**
 *  Solve the eigenvalue problem for a multiplier with the davidson algorithm, davidson guess will be the previously stored solution
 *  
 *  @param multiplier           a given multiplier
 */
void MullikenConstrainedFCI::solveMullikenDavidson(const double multiplier, const double sz_multiplier) {

    if (this->are_solutions_available) {
        this->solveMullikenDavidson(multiplier, eigenvector[0], sz_multiplier);
    } else {
        this->solveMullikenDavidson(multiplier, this->fock_space.hartreeFockExpansion(), sz_multiplier);
    }
}

    
/**
 *  @param index             refers to the index of the number of requested states 
 * 
 *  @return all properties in vector that contains:
 *      energy, population (on the selected basis functions), Sz, lambda (or the multiplier for the Mulliken constraint), lambda_sz (for atomic Sz), entropy
 *      if diatomic we additionally find: A_fragment_energy, A_fragment_self_energy, B_fragment_energy, B_fragment_self_energy and interaction_energy in that order.
 */
void MullikenConstrainedFCI::solveMullikenDense(const double multiplier, const size_t nos = 1, const double sz_multiplier) {
    if (nos < 1 || nos >= fock_space.get_dimension()) {
        throw std::runtime_error("MullikenConstrainedFCI::solveMullikenDense(): number of states should be larger than 0 and smaller than the dimension of the ONV basis.");
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    auto constrained_ham_par = this->usq_hamiltonian.constrain(this->mulliken_operator, multiplier, SpinComponent::ALPHA);
    constrained_ham_par = constrained_ham_par.constrain(this->mulliken_operator, multiplier, SpinComponent::BETA);
    constrained_ham_par = constrained_ham_par.constrain(this->sq_sz_operator, sz_multiplier, SpinComponent::ALPHA);
    constrained_ham_par = constrained_ham_par.constrain(this->sq_sz_operator, -sz_multiplier, SpinComponent::BETA);

    // Dense solver
    const MatrixX<double> H = this->fock_space.evaluateOperatorDense(constrained_ham_par, true);  // the Hamiltonian matrix
    auto dense_environment = GQCP::EigenproblemEnvironment::Dense(H);
    auto dense_diagonalizer = GQCP::EigenproblemSolver::Dense();

    try {
        dense_diagonalizer.perform(dense_environment);
    } catch (const std::exception& e) {
        std::cout << e.what() << "multiplier: " << multiplier;
        return;
    }

    this->parseSolution(dense_environment.eigenpairs(nos), multiplier, sz_multiplier);  // nos: number of requested eigenpairs
    this->are_solutions_available = true;

    auto stop_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = stop_time - start_time;  // in nanoseconds
    this->solve_time = static_cast<double>(elapsed_time.count() / 1e9);  // in seconds
}


std::vector<double> MullikenConstrainedFCI::all_properties(const size_t index) const {
    this->checkAvailableSolutions("all");

    size_t number_of_properties = 5;

    if (this->molecule.numberOfAtoms() == 2) {
        number_of_properties += 6;
    }

    std::vector<double> properties(number_of_properties);

    properties[0] = this->energy[index];
    properties[1] = this->population[index];
    properties[2] = this->sz[index];
    properties[3] = this->lambda[index];
    properties[4] = this->lambda_sz[index];
    properties[5] = this->entropy[index];

    if (molecule.numberOfAtoms() == 2) {
        properties[6] = this->A_fragment_energy[index];
        properties[7] = this->A_fragment_self_energy[index];
        properties[8] = this->B_fragment_energy[index];
        properties[9] = this->B_fragment_self_energy[index];
        properties[10] = this->interaction_energy[index];
    }

    return properties;
}


}  // namespace QCMethod
}  // namespace GQCP


