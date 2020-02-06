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
#define BOOST_TEST_MODULE "DenseDOCISolver"

#include <boost/test/unit_test.hpp>

#include "Basis/transform.hpp"
#include "ONVBasis/ProductONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"


BOOST_AUTO_TEST_CASE ( test_random_rotation_diagonal_dense_fci ) {

    // Check if a random rotation has no effect on the sum of the diagonal elements

    // Create the molecular Hamiltonian in an AO basis
    auto h2o = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2o, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2o);  // in an AO basis
    auto K = sq_hamiltonian.dimension();

    // Create a plain RHF SCF solver and solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2o.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    // Transform the Hamiltonian to the RHF orbital basis
    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());

    GQCP::ProductONVBasis fock_space (K, h2o.numberOfElectrons()/2, h2o.numberOfElectrons()/2);  // dim = 2

    // Create the FCI module
    GQCP::FCI fci (fock_space);

    GQCP::VectorX<double> diagonal1 = fci.calculateDiagonal(sq_hamiltonian);

    // Rotate the hampar using the random unitary matrix
    sq_hamiltonian.randomRotate();

    GQCP::VectorX<double> diagonal2 = fci.calculateDiagonal(sq_hamiltonian);

    BOOST_CHECK(std::abs(diagonal1.sum() - diagonal2.sum()) < 1.0e-10);
}


BOOST_AUTO_TEST_CASE ( FCI_H2_Cristina_dense ) {

    // Cristina's H2 FCI energy/OO-DOCI energy
    double reference_fci_energy = -1.1651486697;

    // Create the molecular Hamiltonian in an AO basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "6-31G**");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis
    auto K = sq_hamiltonian.dimension();

    // Create a plain RHF SCF solver and solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    // Transform the Hamiltonian to the RHF orbital basis
    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());

    GQCP::ProductONVBasis fock_space (K, h2.numberOfElectrons()/2, h2.numberOfElectrons()/2);  // dim = 100

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, sq_hamiltonian);

    // Solve Dense
    GQCP::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2).value();
    double test_fci_energy = fci_energy + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( FCI_H2O_Psi4_GAMESS_dense ) {

    // Psi4 and GAMESS' FCI energy
    double reference_fci_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in an AO basis
    auto h2o = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2o, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2o);  // in an AO basis
    auto K = sq_hamiltonian.dimension();

    // Create a plain RHF SCF solver and solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2o.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    // Transform the Hamiltonian to the RHF orbital basis
    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());

    GQCP::ProductONVBasis fock_space (K, h2o.numberOfElectrons()/2, h2o.numberOfElectrons()/2);  // dim = 441

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, sq_hamiltonian);

    // Solve Dense
    GQCP::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2o).value();
    double test_fci_energy = fci_energy + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}
