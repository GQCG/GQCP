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
#include "QCMethod/Applications/DOCINewtonOrbitalOptimizer.hpp"

#include "Basis/transform.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"
#include "QCMethod/OrbitalOptimization/Localization/ERJacobiLocalizer.hpp"
#include "QCMethod/OrbitalOptimization/Localization/ERNewtonLocalizer.hpp"


namespace GQCP {
namespace QCMethod {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param molecule             the molecule that will be solved for
 *  @param basis_set            the basisset that should be used
 *  @param use_davidson         indicate if one wants to use davidson to solve the eigenvalue problem (opposed to dense)
 *  @param localize             indicate if one wants to localize the orbitals before 
 */
DOCINewtonOrbitalOptimizer::DOCINewtonOrbitalOptimizer(const Molecule& molecule, const std::string& basis_set, const bool use_davidson, const bool localize) :
    molecule (molecule),
    basis_set (basis_set),
    use_davidson (use_davidson),
    localize (localize)
{}


/**
 *  @param xyz_filename         the file that contains the molecule specification (coordinates in angstrom)
 *  @param basis_set            the basisset that should be used
 *  @param use_davidson         indicate if one wants to use davidson to solve the eigenvalue problem (opposed to dense)
 *  @param localize             indicate if one wants to localize the orbitals before 
 */
DOCINewtonOrbitalOptimizer::DOCINewtonOrbitalOptimizer(const std::string& xyz_filename, const std::string& basis_set, const bool use_davidson, const bool localize) :
    DOCINewtonOrbitalOptimizer(Molecule::ReadXYZ(xyz_filename), basis_set, use_davidson, localize)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the dense eigenvalue problem for the molecular Hamiltonian in the full ONV basis
 */
void DOCINewtonOrbitalOptimizer::solve() {

    // Construct the molecular Hamiltonian in the RHF basis
    RSpinorBasis<double, GTOShell> spinor_basis (this->molecule, this->basis_set);
    const auto N_P = molecule.numberOfElectrons()/2;
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, this->molecule);  // in AO basis
    const size_t K = sq_hamiltonian.dimension();

    // Do an RHF calculation and transform the spinor basis to the RHF basis
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(this->molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, rhf_scf_solver, rhf_environment).groundStateParameters();

    basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());

    auto hessian_modifier = std::make_shared<IterativeIdentitiesHessianModifier>();
    if (localize) {

        // Newton to get to the first local minimum
        ERNewtonLocalizer first_newton_localizer (N_P, hessian_modifier);
        first_newton_localizer.optimize(spinor_basis, sq_hamiltonian);


        // Check if Jacobi finds another minimum
        ERJacobiLocalizer jacobi_localizer (N_P);
        auto optimal_jacobi_with_scalar = jacobi_localizer.calculateOptimalJacobiParameters(sq_hamiltonian);
        if (optimal_jacobi_with_scalar.second > 0) {  // if a Jacobi rotation can find an increase, do it
            const TransformationMatrix<double> U = TransformationMatrix<double>::FromJacobi(optimal_jacobi_with_scalar.first, sq_hamiltonian.dimension());
            basisRotate(spinor_basis, sq_hamiltonian, U);
        }


        // Newton to get to the next local minimum
        ERNewtonLocalizer second_newton_localizer (N_P, hessian_modifier);
        first_newton_localizer.optimize(spinor_basis, sq_hamiltonian);
    }

     // Do the DOCI orbital optimization
    SpinUnresolvedONVBasis fock_space (K, N_P);
    DOCI doci (fock_space);

    std::shared_ptr<BaseSolverOptions> solver_options;
    if (use_davidson) {
        solver_options = std::make_shared<DavidsonSolverOptions>(fock_space.HartreeFockExpansion());
    } else {
        solver_options = std::make_shared<DenseSolverOptions>();
    }

    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, *solver_options, hessian_modifier);
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);
    double OO_DOCI_electronic_energy = orbital_optimizer.get_eigenpair().get_eigenvalue();

    this->is_solved = true;
    double internuclear_repulsion_energy = Operator::NuclearRepulsion(molecule).value();
    this->energy_solution = OO_DOCI_electronic_energy + internuclear_repulsion_energy;
    this->T_total = spinor_basis.coefficientMatrix();
}


/**
 *  @return the newton orbital optimized ground state DOCI energy
 */  
double DOCINewtonOrbitalOptimizer::energy() const {

    if (this->is_solved) {
        return this->energy_solution;
    } else {
        throw std::runtime_error("DOCINewtonOrbitalOptimizer::energy(): You are trying to get energy but the method hasn't been solved yet.");
    }
}


/**
 *  @return the total transformation matrix to the OO-DOCI orbitals
 */
const TransformationMatrix<double>& DOCINewtonOrbitalOptimizer::transformationMatrix() const {

    if (this->is_solved) {
        return this->T_total;
    } else {
        throw std::runtime_error("DOCINewtonOrbitalOptimizer::transformationMatrix(): You are trying to get the total transformation matrix but the method hasn't been solved yet.");
    }
}


}  // namespace QCMethod
}  // namespace GQCP
