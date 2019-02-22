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
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "CISolver/CISolver.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"



BOOST_AUTO_TEST_CASE ( test_random_rotation_diagonal_dense_fci ) {

    // Check if a random rotation has no effect on the sum of the diagonal elements

    // Create the molecular Hamiltonian parameters in an AO basis
    auto h2o = GQCP::Molecule::Readxyz("data/h2o.xyz");
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(h2o, "STO-3G");
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 2

    // Create the FCI module
    GQCP::FCI fci (fock_space);

    Eigen::VectorXd diagonal1 = fci.calculateDiagonal(mol_ham_par);

    // Get a random unitary matrix by diagonalizing a random symmetric matrix
    Eigen::MatrixXd A_random = Eigen::MatrixXd::Random(K, K);
    Eigen::MatrixXd A_symmetric = A_random + A_random.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> unitary_solver (A_symmetric);
    Eigen::MatrixXd U_random = unitary_solver.eigenvectors();

    // Rotate the hampar using the random unitary matrix
    mol_ham_par.rotate(U_random);

    Eigen::VectorXd diagonal2 = fci.calculateDiagonal(mol_ham_par);

    BOOST_CHECK(std::abs(diagonal1.sum() - diagonal2.sum()) < 1.0e-10);
}


BOOST_AUTO_TEST_CASE ( FCI_H2_Cristina_dense ) {

    // Cristina's H2 FCI energy/OO-DOCI energy
    double reference_fci_energy = -1.1651486697;

    // Create the molecular Hamiltonian parameters in an AO basis
    auto h2 = GQCP::Molecule::Readxyz("data/h2_cristina.xyz");
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(h2, "6-31g**");
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2.get_N()/2, h2.get_N()/2);  // dim = 100

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, mol_ham_par);

    // Solve Dense
    GQCP::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci_energy + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( FCI_H2O_Psi4_GAMESS_dense ) {

    // Psi4 and GAMESS' FCI energy
    double reference_fci_energy = -75.0129803939602;

    // Create the molecular Hamiltonian parameters in an AO basis
    auto h2o = GQCP::Molecule::Readxyz("data/h2o_Psi4_GAMESS.xyz");
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(h2o, "STO-3G");
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 441

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, mol_ham_par);

    // Solve Dense
    GQCP::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = h2o.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci_energy + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}

/*
BOOST_AUTO_TEST_CASE ( FCI_He_Cristina_dense ) {

    // Cristina's He FCI energy
    double reference_fci_energy = -2.902533599;

    // Create a Molecule and an AOBasis
    GQCP::Molecule he ("data/h2o_Psi4_GAMESS.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(he, "aug-cc-pVQZ");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, he);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, he.get_N()/2, he.get_N()/2);  // dim = 2116

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, mol_ham_par);

    // Solve Dense
    GQCP::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = he.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci_energy + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}
 */
