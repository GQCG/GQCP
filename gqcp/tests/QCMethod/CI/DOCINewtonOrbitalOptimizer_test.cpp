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
#define BOOST_TEST_MODULE "DOCI_orbital_optimization_test"

#include <boost/test/unit_test.hpp>

#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"

#include "Basis/transform.hpp"
#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/FCIRDMBuilder.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"


/**
 *  Check if OO-DOCI (dense) matches FCI for a two-electron system.
 *  The system of interested is H2//STO-3G, with reference results obtained from Christina at Ayer's lab.
 */
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_sto_3g ) {

    const double reference_fci_energy = -1.13726333769813;

    // Prepare the molecular Hamiltonian in the canonical RHF basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;
    const auto internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(molecule).value();  // 0.713176780299327

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Do the DOCI orbital optimization.
    const GQCP::SeniorityZeroONVBasis onv_basis (K, N_P);

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    using EigenproblemSolver = decltype(solver);

    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    GQCP::DOCINewtonOrbitalOptimizer<EigenproblemSolver> orbital_optimizer (onv_basis, solver, environment, hessian_modifier);
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);

    const auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();


    // Check if the OO-DOCI energy is equal to the FCI energy.
    const double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}


/**
 *  Check if OO-DOCI (Davidson) matches FCI for a two-electron system.
 *  The system of interested is H2//6-31G**, with reference results obtained from Christina at Ayer's lab.
 */
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31gxx_Davidson ) {

    const double reference_fci_energy = -1.16514875501195;

    // Prepare the molecular Hamiltonian in the canonical RHF basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;
    const auto internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(molecule).value();  // 0.713176780299327

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "6-31G**");
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Do the DOCI orbital optimization.
    const GQCP::SeniorityZeroONVBasis onv_basis (K, N_P);

    const auto initial_guess = onv_basis.hartreeFockExpansion();
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, initial_guess);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    using EigenproblemSolver = decltype(solver);





    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    GQCP::DOCINewtonOrbitalOptimizer<EigenproblemSolver> orbital_optimizer (onv_basis, solver, environment, hessian_modifier);
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);

    const auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();


    // Check if the OO-DOCI energy is equal to the FCI energy.
    const double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    std::cout << "OO-DOCI-energy " << OO_DOCI_energy << std::endl;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}
