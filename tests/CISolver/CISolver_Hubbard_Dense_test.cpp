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


#include "CISolver/CISolver.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( test_random_rotation_diagonal_dense_Hubbard ) {

    // Check if a random rotation has no effect on the sum of the diagonal elements


    // Create the Hamiltonian parameters for the triagonal of a Hubbard lattice.
    Eigen::VectorXd triagonal_test(6);
    triagonal_test << 1, 2, 3, 4, 5, 6;

    size_t N = 2;
    auto mol_ham_par = GQCP::hubbardTriagonalLattice(triagonal_test);
    auto K = mol_ham_par.get_K();


    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 9

    // Create the Hubbard module
    GQCP::Hubbard hubbard (fock_space);

    Eigen::VectorXd diagonal1 = hubbard.calculateDiagonal(mol_ham_par);

    // Get a random unitary matrix by diagonalizing a random symmetric matrix
    Eigen::MatrixXd A_random = Eigen::MatrixXd::Random(K, K);
    Eigen::MatrixXd A_symmetric = A_random + A_random.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> unitary_solver (A_symmetric);
    Eigen::MatrixXd U_random = unitary_solver.eigenvectors();

    // Rotate the hampar using the random unitary matrix
    mol_ham_par.rotate(U_random);

    Eigen::VectorXd diagonal2 = hubbard.calculateDiagonal(mol_ham_par);

    BOOST_CHECK(std::abs(diagonal1.sum() - diagonal2.sum()) < 1.0e-10);
}

BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_dense ) {

    // Check if FCI and Hubbard produce the same results

    // Create the Hamiltonian parameters for the triagonal of a Hubbard lattice.
    Eigen::VectorXd triagonal_test = Eigen::VectorXd::Random(10);

    size_t N = 2;
    auto mol_ham_par = GQCP::hubbardTriagonalLattice(triagonal_test);
    auto K = mol_ham_par.get_K();


    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 36

    // Create the Hubbard module
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    GQCP::CISolver solver1 (hubbard, mol_ham_par);
    GQCP::CISolver solver2 (fci, mol_ham_par);

    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    solver1.solve(dense_solver_options);
    solver2.solve(dense_solver_options);

    auto fci_energy = solver2.get_eigenpair().get_eigenvalue();
    auto hubbard_energy = solver1.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_energy - (hubbard_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_dense_large ) {

    // Check if FCI and Hubbard produce the same results

    // Create the Hamiltonian parameters for the triagonal of a Hubbard lattice.
    Eigen::VectorXd triagonal_test = Eigen::VectorXd::Random(21);

    size_t N = 3;
    auto mol_ham_par = GQCP::hubbardTriagonalLattice(triagonal_test);
    auto K = mol_ham_par.get_K();


    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 400

    // Create the Hubbard module
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    GQCP::CISolver solver1 (hubbard, mol_ham_par);
    GQCP::CISolver solver2 (fci, mol_ham_par);

    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    solver1.solve(dense_solver_options);
    solver2.solve(dense_solver_options);

    auto fci_energy = solver2.get_eigenpair().get_eigenvalue();
    auto hubbard_energy = solver1.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_energy - (hubbard_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( four_site_chain_ward ) {

    size_t K = 4;

    // Create the adjacency matrix
    Eigen::MatrixXd A (K, K);
    A << 0, 1, 0, 0,
         1, 0, 1, 0,
         0, 1, 0, 1,
         0, 0, 1, 0;

    double t = 1.0;
    std::vector<double> U_list {0.0, 1.0, 1.5, 3.5, 6.5, 9, 10};
    size_t N = 2;

    // Set Ward's results
    std::vector<double> E_list {-4.472135955, -3.57536562, -3.202271824, -2.135871608, -1.338959715, -1.004379799, -0.9114974686};

    // Specify a tolerance
    double tol = 1.0e-08;

    GQCP::ProductFockSpace fock_space (K, N, N);

    // For every U-value, create the Hubbard instance, solve it, and check the energy with Ward's results (https://github.com/wpoely86/Hubbard-GPU)
    for (size_t i = 0; i < 7; i++) {

        GQCP::Hubbard hubbard (fock_space);
        Eigen::VectorXd triagonal = GQCP::genrateUpperTriagonal(A, t, U_list[i]);
        auto ham_par = GQCP::hubbardTriagonalLattice(triagonal);

        GQCP::CISolver solver1 (hubbard, ham_par);
        numopt::eigenproblem::DenseSolverOptions dense_solver_options;
        solver1.solve(dense_solver_options);

        BOOST_CHECK(std::abs(solver1.get_eigenpair().get_eigenvalue() - E_list[i]) < tol);
    }
}


BOOST_AUTO_TEST_CASE ( six_site_ring_ward ) {

    size_t K = 6;

    // Create the adjacency matrix
    Eigen::MatrixXd A (K, K);
    A << 0, 1, 0, 0, 0, 1,
         1, 0, 1, 0, 0, 0,
         0, 1, 0, 1, 0, 0,
         0, 0, 1, 0, 1, 0,
         0, 0, 0, 1, 0, 1,
         1, 0, 0, 0, 1, 0;

    double t = 1.0;

    // Set the HubbardClass parameters
    size_t N = 6;

    std::vector<double> U_list {0.0, 1.0, 1.5, 3.5, 6.5, 9, 10};

    // Set Ward's results
    std::vector<double> E_list {-8, -6.601158293, -5.978815789, -4.025796251, -2.469458295, -1.836926909, -1.664362733};

    // Specify a tolerance
    double tol = 1.0e-08;

    GQCP::ProductFockSpace fock_space (K, N, N);

    // For every U-value, create the Hubbard instance, solve it, and check the energy with Ward's results (https://github.com/wpoely86/Hubbard-GPU)
    for (size_t i = 0; i < 7; i++) {

        GQCP::Hubbard hubbard (fock_space);
        Eigen::VectorXd triagonal = GQCP::genrateUpperTriagonal(A, t, U_list[i]);
        auto ham_par = GQCP::hubbardTriagonalLattice(triagonal);

        GQCP::CISolver solver1 (hubbard, ham_par);
        numopt::eigenproblem::DenseSolverOptions dense_solver_options;
        solver1.solve(dense_solver_options);

        BOOST_CHECK(std::abs(solver1.get_eigenpair().get_eigenvalue() - E_list[i]) < tol);
    }
}

