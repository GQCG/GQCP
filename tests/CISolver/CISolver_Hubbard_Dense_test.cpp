// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#include "HamiltonianBuilder/Hubbard.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_dense ) {

    // Check if FCI and Hubbard produce the same results for Hubbard Hamiltonian parameters

    // Create the Hamiltonian parameters for a random Hubbard hopping matrix
    size_t K = 4;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto mol_ham_par = GQCP::HamiltonianParameters::Hubbard(H);


    // Create the Hubbard and FCI modules
    size_t N = 2;
    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 36
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);


    // Solve via dense
    GQCP::CISolver hubbard_solver (hubbard, mol_ham_par);
    GQCP::CISolver fci_solver (fci, mol_ham_par);

    GQCP::DenseSolverOptions dense_solver_options;
    hubbard_solver.solve(dense_solver_options);
    fci_solver.solve(dense_solver_options);

    auto fci_energy = fci_solver.get_eigenpair().get_eigenvalue();
    auto hubbard_energy = hubbard_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_energy - (hubbard_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_dense_large ) {

    // Check if FCI and Hubbard produce the same results for Hubbard Hamiltonian parameters

    // Create the Hamiltonian parameters for a random Hubbard hopping matrix
    size_t K = 6;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto mol_ham_par = GQCP::HamiltonianParameters::Hubbard(H);


    // Create the Hubbard and FCI modules
    size_t N = 3;
    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 36
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    GQCP::CISolver solver1 (hubbard, mol_ham_par);
    GQCP::CISolver solver2 (fci, mol_ham_par);

    GQCP::DenseSolverOptions dense_solver_options;
    solver1.solve(dense_solver_options);
    solver2.solve(dense_solver_options);

    auto fci_energy = solver2.get_eigenpair().get_eigenvalue();
    auto hubbard_energy = solver1.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_energy - (hubbard_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( four_site_chain_ward ) {

    // Create the adjacency matrix for a four-site chain
    size_t K = 4;
    size_t N = 2;  // = N_alpha = N_beta: half-filling
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(K, K);
    A << 0, 1, 0, 0,
         1, 0, 1, 0,
         0, 1, 0, 1,
         0, 0, 1, 0;

    // Set Ward's results (https://github.com/wpoely86/Hubbard-GPU) from ground state energies
    double t = 1.0;
    std::vector<double> U_list {0.0, 1.0, 1.5, 3.5, 6.5, 9, 10};
    std::vector<double> E_list {-4.472135955, -3.57536562, -3.202271824, -2.135871608, -1.338959715, -1.004379799, -0.9114974686};
    double tol = 1.0e-08;


    GQCP::ProductFockSpace fock_space (K, N, N);
    for (size_t i = 0; i < 7; i++) {

        // Create the Hamiltonian parameters for the Hubbard model
        GQCP::HoppingMatrix H (A, t, U_list[i]);
        auto ham_par = GQCP::HamiltonianParameters::Hubbard(H);


        // Solve the dense eigenvalue problem
        GQCP::Hubbard hubbard (fock_space);
        GQCP::CISolver solver (hubbard, ham_par);
        GQCP::DenseSolverOptions dense_solver_options;
        solver.solve(dense_solver_options);

        BOOST_CHECK(std::abs(solver.get_eigenpair().get_eigenvalue() - E_list[i]) < tol);
    }
}


BOOST_AUTO_TEST_CASE ( six_site_ring_ward ) {

    // Create the adjacency matrix for a six-site ring
    size_t K = 6;
    size_t N = 3;  // = N_alpha = N_beta: half-filling
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(K, K);
    A << 0, 1, 0, 0, 0, 1,
         1, 0, 1, 0, 0, 0,
         0, 1, 0, 1, 0, 0,
         0, 0, 1, 0, 1, 0,
         0, 0, 0, 1, 0, 1,
         1, 0, 0, 0, 1, 0;


    // Set Ward's results (https://github.com/wpoely86/Hubbard-GPU) from ground state energies
    double t = 1.0;
    std::vector<double> U_list {0.0, 1.0, 1.5, 3.5, 6.5, 9, 10};
    std::vector<double> E_list {-8, -6.601158293, -5.978815789, -4.025796251, -2.469458295, -1.836926909, -1.664362733};
    double tol = 1.0e-08;


    GQCP::ProductFockSpace fock_space (K, N, N);
    for (size_t i = 0; i < 7; i++) {

        // Create the Hamiltonian parameters for the Hubbard model
        GQCP::HoppingMatrix H (A, t, U_list[i]);
        auto ham_par = GQCP::HamiltonianParameters::Hubbard(H);


        // Solve the dense eigenvalue problem
        GQCP::Hubbard hubbard (fock_space);
        GQCP::CISolver solver (hubbard, ham_par);
        GQCP::DenseSolverOptions dense_solver_options;
        solver.solve(dense_solver_options);

        BOOST_CHECK(std::abs(solver.get_eigenpair().get_eigenvalue() - E_list[i]) < tol);
    }
}

