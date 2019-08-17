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
#include "QCMethod/MullikenConstrainedFCI.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"
#include "Properties/expectation_values.hpp"

namespace GQCP {
namespace QCMethod {

/*
 *  PRIVATE METHODS
 */

/**
 *  Store the solutions from a solve
 *  
 *  @param eigenpairs           the eigenpairs from the CI solver
 */ 
void MullikenConstrainedFCI::parseSolution(const std::vector<Eigenpair>& eigenpairs, double multiplier) {
    if (this->energy.size() != eigenpairs.size()) {

        this->energy = std::vector<double>(eigenpairs.size());
        this->population = std::vector<double>(eigenpairs.size());
        this->lambda = std::vector<double>(eigenpairs.size());
        this->entropy = std::vector<double>(eigenpairs.size());

        if (molecule.numberOfAtoms() == 2) {
            this->A_fragment_energy = std::vector<double>(eigenpairs.size());
            this->A_fragment_self_energy = std::vector<double>(eigenpairs.size());
            this->B_fragment_energy = std::vector<double>(eigenpairs.size());
            this->B_fragment_self_energy = std::vector<double>(eigenpairs.size());
            this->interaction_energy = std::vector<double>(eigenpairs.size());
        }

        this->eigenvector = std::vector<VectorX<double>>(eigenpairs.size());
    }

    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(molecule).value();

    for (size_t i = 0; i < eigenpairs.size(); i++) {

        const auto& pair = eigenpairs[i];
        const auto& fci_coefficients = pair.get_eigenvector();
        double fci_energy = pair.get_eigenvalue();
        rdm_calculator.set_coefficients(fci_coefficients);
        OneRDM<double> D = rdm_calculator.calculate1RDMs().one_rdm;
        TwoRDM<double> d = rdm_calculator.calculate2RDMs().two_rdm;

        double population = calculateExpectationValue(mulliken_operator, D);
        WaveFunction wavefunction (fock_space, fci_coefficients);

        this->energy[i] = pair.get_eigenvalue() + internuclear_repulsion_energy + multiplier * population;
        this->population[i] = population;
        this->lambda[i] = multiplier;
        this->entropy[i] = wavefunction.calculateShannonEntropy();

        if (molecule.numberOfAtoms() == 2) {
            this->A_fragment_energy[i] = calculateExpectationValue(adp.get_atomic_parameters()[0], D, d);
            this->A_fragment_self_energy[i] = calculateExpectationValue(adp.get_net_atomic_parameters()[0], D, d);
            this->B_fragment_energy[i] = calculateExpectationValue(adp.get_atomic_parameters()[1], D, d);
            this->B_fragment_self_energy[i] = calculateExpectationValue(adp.get_net_atomic_parameters()[1], D, d);
            this->interaction_energy[i] = calculateExpectationValue(adp.get_interaction_parameters()[0], D, d);
        }

        this->eigenvector[i] = fci_coefficients;
    }
}


/**
 *  Throws and error if no solution is available
 *  
 *  @param function_name            name of the function that should throw the error
 */
void MullikenConstrainedFCI::checkAvailableSolutions(const std::string& function_name) const {
    if (!are_solutions_available) {
        throw std::runtime_error("MullikenConstrainedFCI::" + function_name + "(): The method hasn't been solved yet");
    }
}

/**
 *  Throws and error if the molecule is not diatomic
 *  
 *  @param function_name            name of the function that should throw the error
 */
void MullikenConstrainedFCI::checkDiatomicMolecule(const std::string& function_name) const {
    if (molecule.numberOfAtoms() != 2) {
        throw std::runtime_error("MullikenConstrainedFCI::" + function_name + "(): This property only available for diatomic molecules");
    }
}

/*
 * CONSTRUCTORS
 */

/**
 *  @param molecule                 the molecule that will be solved for
 *  @param basis_set                the basisset that should be used
 *  @param basis_targets            the targeted basis functions for the constraint
 *  @param multipliers              the set of multipliers for the constraint
 */
MullikenConstrainedFCI::MullikenConstrainedFCI(const Molecule& molecule, const std::string& basis_set, const std::vector<size_t>& basis_targets, size_t frozencores) : 
        basis_targets (basis_targets),
        molecule (molecule),
        ham_par (HamiltonianParameters<double>::Molecular(molecule, basis_set)),
        basis_set (basis_set)
{
    if (!(molecule.numberOfElectrons() % 2)) {
        throw std::runtime_error("MullikenConstrainedFCI::MullikenConstrainedFCI(): This module is not available for an odd number of electrons");
    }

    auto K = this->ham_par.get_K();
    auto N_P = molecule.numberOfElectrons()/2;
    try {
        GQCP::DIISRHFSCFSolver diis_scf_solver (this->ham_par, molecule, 6, 6, 1e-12, 500);
        auto rhf_solution = diis_scf_solver.get_solution();
        this->ham_par.basisTransform(rhf_solution.get_C());
    } catch (const std::exception& e) {
        
        if (molecule.numberOfAtoms() == 2) {
            try {
                const std::vector<Nucleus>& atoms = molecule.nuclearFramework().nucleiAsVector();
                int charge = - molecule.numberOfElectrons() + molecule.nuclearFramework().totalNucleicCharge();
                GQCP::Molecule mol_fraction1(std::vector<GQCP::Nucleus>{atoms[0]}, charge);
                GQCP::Molecule mol_fraction2(std::vector<GQCP::Nucleus>{atoms[1]}, 0);

                auto ham_par1 = GQCP::HamiltonianParameters<double>::Molecular(mol_fraction1, basis_set);
                auto ham_par2 = GQCP::HamiltonianParameters<double>::Molecular(mol_fraction2, basis_set);

                // Perform DIIS RHF for individual fractions
                GQCP::DIISRHFSCFSolver diis_scf_solver1 (ham_par1, mol_fraction1, 6, 6, 1e-12, 500);
                GQCP::DIISRHFSCFSolver diis_scf_solver2 (ham_par2, mol_fraction2, 6, 6, 1e-12, 500);
                diis_scf_solver1.solve();
                diis_scf_solver2.solve();
                auto rhf1 = diis_scf_solver1.get_solution();
                auto rhf2 = diis_scf_solver2.get_solution();

                // Retrieve transformation from the solutions and transform the Hamiltonian parameters
                size_t K1 = ham_par1.get_K();
                size_t K2 = ham_par2.get_K();

                GQCP::SquareMatrix<double> T = Eigen::MatrixXd::Zero(K, K);
                T.topLeftCorner(K1, K1) += rhf1.get_C();
                T.bottomRightCorner(K2, K2) += rhf2.get_C();
                this->ham_par.basisTransform(T);

                try {
                    GQCP::DIISRHFSCFSolver diis_scf_solver (this->ham_par, molecule, 6, 6, 1e-12, 500);
                    diis_scf_solver.solve();
                    auto rhf = diis_scf_solver.get_solution();
                    this->ham_par.basisTransform(rhf.get_C());

                } catch (const std::exception& e) {
                    std::cout << "Lodwin Orthonormalized" << std::endl;
                    this->ham_par.LowdinOrthonormalize();
                }
            } catch (const std::exception& e) {
                std::cout << "Lodwin Orthonormalized" << std::endl;
                this->ham_par.LowdinOrthonormalize();
            }

        } else {
            std::cout << "Lodwin Orthonormalized" << std::endl;
            this->ham_par.LowdinOrthonormalize();
        }

    }

    this->fock_space = FrozenProductFockSpace(K, N_P, N_P, frozencores);
    this->fci = FrozenCoreFCI(fock_space);
    this->mulliken_operator = ham_par.calculateMullikenOperator(basis_targets);

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
 */
void MullikenConstrainedFCI::solveMullikenDavidson(const double multiplier, const VectorX<double>& guess) {

    auto start_time = std::chrono::high_resolution_clock::now();

    const auto constrained_ham_par = this->ham_par.constrain(mulliken_operator, multiplier);
    CISolver ci_solver (fci, constrained_ham_par);
    DavidsonSolverOptions solver_options (fock_space.HartreeFockExpansion());
    try {
        ci_solver.solve(solver_options);
    } catch (const std::exception& e) {
        std::cout << e.what() << "multiplier: " << multiplier;
        return;
    }

    this->parseSolution(ci_solver.get_eigenpairs(), multiplier);
    this->are_solutions_available = true;

    auto stop_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = stop_time- start_time;  // in nanoseconds
    this->solve_time = static_cast<double>(elapsed_time.count() / 1e9);  // in seconds
}

/**
 *  Solve the eigenvalue problem for a multiplier with the davidson algorithm, davidson guess will be the previously stored solution
 *  
 *  @param multiplier           a given multiplier          
 */
void MullikenConstrainedFCI::solveMullikenDavidson(const double multiplier) {

    if (this->are_solutions_available) {
        this->solveMullikenDavidson(multiplier, eigenvector[0]);
    } else {
        this->solveMullikenDavidson(multiplier, this->fock_space.HartreeFockExpansion());
    }
}

/**
 *  Solve the eigenvalue problem for a the next multiplier dense
 * 
 *  @param multiplier     
 *  @param nos                  the number of eigenpairs or "states" that should be stored for each multiplier
 */
void MullikenConstrainedFCI::solveMullikenDense(const double multiplier, const size_t nos = 1) {
    if (nos < 1 || nos >= fock_space.get_dimension()) {
        throw std::runtime_error("MullikenConstrainedFCI::solveMullikenDense(): number of states should be larger than 0 and smaller than the dimension of the Fock space");
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    const auto constrained_ham_par = this->ham_par.constrain(mulliken_operator, multiplier);
    CISolver ci_solver (fci, constrained_ham_par);
    DenseSolverOptions solver_options;
    solver_options.number_of_requested_eigenpairs = nos;

    try {
        ci_solver.solve(solver_options);
    } catch (const std::exception& e) {
        std::cout << e.what() << "multiplier: " << multiplier;
        return;
    }

    this->parseSolution(ci_solver.get_eigenpairs(), multiplier);
    this->are_solutions_available = true;

    auto stop_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = stop_time- start_time;  // in nanoseconds
    this->solve_time = static_cast<double>(elapsed_time.count() / 1e9);  // in seconds
}


std::vector<double> MullikenConstrainedFCI::all(size_t index) const {
    this->checkAvailableSolutions("all");

    size_t number_of_properties = 4;

    if (molecule.numberOfAtoms() == 2) {
        number_of_properties += 5;
    }

    std::vector<double> properties(number_of_properties);

    properties[0] = this->energy[index];
    properties[1] = this->population[index];
    properties[2] = this->lambda[index];
    properties[3] = this->entropy[index];

    if (molecule.numberOfAtoms() == 2) {
        properties[4] = this->A_fragment_energy[index];
        properties[5] = this->A_fragment_self_energy[index];
        properties[6] = this->B_fragment_energy[index];
        properties[7] = this->B_fragment_self_energy[index];
        properties[8] = this->interaction_energy[index];
    }

    return properties;
}


}  // namespace QCMethod
}  // namespace GQCP


