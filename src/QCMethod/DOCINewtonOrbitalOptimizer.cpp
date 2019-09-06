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
#include "QCMethod/DOCINewtonOrbitalOptimizer.hpp"

#include "CISolver/CISolver.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Mathematical/Optimization/IterativeIdentitiesHessianModifier.hpp"
#include "OrbitalOptimization/DOCINewtonOrbitalOptimizer.hpp"
#include "OrbitalOptimization/Localization/ERJacobiLocalizer.hpp"
#include "OrbitalOptimization/Localization/ERNewtonLocalizer.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"


namespace GQCP {
namespace QCMethod {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param xyz_filename         the file that contains the molecule specification (coordinates in angstrom)
 *  @param basis_set            the basisset that should be used
 *  @param use_davidson         indicate if one wants to use davidson to solve the eigenvalue problem (opposed to dense)
 *  @param localize             indicate if one wants to localize the orbitals before 
 */
DOCINewtonOrbitalOptimizer::DOCINewtonOrbitalOptimizer(const std::string& xyz_filename, const std::string& basis_set, const bool use_davidson, const bool localize) :
    xyz_filename (xyz_filename),
    basis_set (basis_set),
    use_davidson (use_davidson),
    localize (localize)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the dense eigenvalue problem for the molecular Hamiltonian in the full Fock space
 */
void DOCINewtonOrbitalOptimizer::solve() {

    // Construct the molecular Hamiltonian parameters
    auto molecule = Molecule::ReadXYZ(this->xyz_filename);
    auto N_P = molecule.numberOfElectrons()/2;
    auto ao_mol_ham_par = GQCP::SQHamiltonian<double>::Molecular(molecule, basis_set);
    size_t K = ao_mol_ham_par.get_K();


    GQCP::DIISRHFSCFSolver diis_scf_solver (ao_mol_ham_par, molecule);
    diis_scf_solver.solve();
    auto rhf = diis_scf_solver.get_solution();

    auto mol_ham_par = GQCP::SQHamiltonian<double>(ao_mol_ham_par, rhf.get_C());

    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    if (localize) {

        // Newton to get to the first local minimum
        GQCP::ERNewtonLocalizer first_newton_localizer (N_P, hessian_modifier);
        first_newton_localizer.optimize(mol_ham_par);


        // Check if Jacobi finds another minimum
        GQCP::ERJacobiLocalizer jacobi_localizer (N_P);
        auto optimal_jacobi_with_scalar = jacobi_localizer.calculateOptimalJacobiParameters(mol_ham_par);
        if (optimal_jacobi_with_scalar.second > 0) {  // if a Jacobi rotation can find an increase, do it
            const auto U = GQCP::TransformationMatrix<double>::FromJacobi(optimal_jacobi_with_scalar.first, mol_ham_par.get_K());
            mol_ham_par.rotate(U);
        }


        // Newton to get to the next local minimum
        GQCP::ERNewtonLocalizer second_newton_localizer (N_P, hessian_modifier);
        first_newton_localizer.optimize(mol_ham_par);
    }

     // Do the DOCI orbital optimization
    GQCP::FockSpace fock_space (K, N_P);
    GQCP::DOCI doci (fock_space);

    std::shared_ptr<GQCP::BaseSolverOptions> solver_options;
    if (use_davidson) {
        solver_options = std::make_shared<GQCP::DavidsonSolverOptions>(fock_space.HartreeFockExpansion());
    } else {
        solver_options = std::make_shared<GQCP::DenseSolverOptions>();
    }

    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, *solver_options, hessian_modifier);
    orbital_optimizer.optimize(mol_ham_par);
    double OO_DOCI_electronic_energy = orbital_optimizer.get_eigenpair().get_eigenvalue();

    this->is_solved = true;
    double internuclear_repulsion_energy = Operator::NuclearRepulsion(molecule).value();
    this->energy_solution = OO_DOCI_electronic_energy + internuclear_repulsion_energy;
    this->T_total = mol_ham_par.get_T_total();  
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
const GQCP::TransformationMatrix<double>& DOCINewtonOrbitalOptimizer::transformationMatrix() const {

    if (this->is_solved) {
        return this->T_total;
    } else {
        throw std::runtime_error("DOCINewtonOrbitalOptimizer::transformationMatrix(): You are trying to get the total transformation matrix but the method hasn't been solved yet.");
    }
}


}  // namespace QCMethod
}  // namespace GQCP
