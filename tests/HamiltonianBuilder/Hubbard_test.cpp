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
#define BOOST_TEST_MODULE "Hubbard"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "HamiltonianBuilder/Hubbard.hpp"

#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"



BOOST_AUTO_TEST_CASE ( Hubbard_constructor ) {

    // Check if a correct constructor works
    GQCP::ProductFockSpace fock_space (15, 3, 3);
    BOOST_CHECK_NO_THROW(GQCP::Hubbard Hubbard (fock_space));
}


BOOST_AUTO_TEST_CASE ( Hubbard_public_methods ) {

    size_t K = 7;
    auto random_hamiltonian_parameters = GQCP::HamiltonianParameters<double>::Random(K);


    // Create a compatible Fock space
    GQCP::ProductFockSpace fock_space (K, 3, 3);
    GQCP::Hubbard random_Hubbard (fock_space);
    Eigen::VectorXd x = random_Hubbard.calculateDiagonal(random_hamiltonian_parameters);
    BOOST_CHECK_NO_THROW(random_Hubbard.constructHamiltonian(random_hamiltonian_parameters));
    BOOST_CHECK_NO_THROW(random_Hubbard.matrixVectorProduct(random_hamiltonian_parameters, x, x));


    // Create an incompatible Fock space
    GQCP::ProductFockSpace fock_space_invalid (K+1, 3, 3);

    // Create Hubbard module
    GQCP::Hubbard random_Hubbard_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_Hubbard_invalid.constructHamiltonian(random_hamiltonian_parameters), std::invalid_argument);
    BOOST_CHECK_THROW(random_Hubbard_invalid.matrixVectorProduct(random_hamiltonian_parameters, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI ) {

    // Create Hubbard Hamiltonian parameters
    size_t K = 4;
    size_t N = 2;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Hubbard(H);


    // Check if the Hamiltonian matrix is equal for FCI and Hubbard
    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 36
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    Eigen::MatrixXd hubbard_ham = hubbard.constructHamiltonian(mol_ham_par);
    Eigen::MatrixXd fci_ham = fci.constructHamiltonian(mol_ham_par);

    BOOST_CHECK(hubbard_ham.isApprox(fci_ham));
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_large ) {

    // Create Hubbard Hamiltonian parameters
    size_t K = 6;
    size_t N = 3;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Hubbard(H);


    // Check if the Hamiltonian matrix is equal for FCI and Hubbard
    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 400
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    Eigen::MatrixXd hubbard_ham = hubbard.constructHamiltonian(mol_ham_par);
    Eigen::MatrixXd fci_ham = fci.constructHamiltonian(mol_ham_par);
    BOOST_CHECK(hubbard_ham.isApprox(fci_ham));
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_matvec ) {

    // Create Hubbard Hamiltonian parameters
    size_t K = 4;
    size_t N = 2;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Hubbard(H);


    // Check if FCI and Hubbard have the same matvec for Hubbard Hamiltonian parameters
    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 36

    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    Eigen::VectorXd hubbard_diagonal = hubbard.calculateDiagonal(mol_ham_par);
    Eigen::VectorXd fci_diagonal = fci.calculateDiagonal(mol_ham_par);
    BOOST_CHECK(hubbard_diagonal.isApprox(fci_diagonal));

    Eigen::VectorXd hubbard_matvec = hubbard.matrixVectorProduct(mol_ham_par, hubbard_diagonal, hubbard_diagonal);
    Eigen::VectorXd fci_matvec = fci.matrixVectorProduct(mol_ham_par, fci_diagonal, fci_diagonal);
    BOOST_CHECK(hubbard_matvec.isApprox(fci_matvec));
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_large_matvec ) {

    // Create Hubbard Hamiltonian parameters
    size_t K = 6;
    size_t N = 3;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Hubbard(H);


    // Check if FCI and Hubbard have the same matvec for Hubbard Hamiltonian parameters
    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 400

    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    Eigen::VectorXd hubbard_diagonal = hubbard.calculateDiagonal(mol_ham_par);
    Eigen::VectorXd fci_diagonal = fci.calculateDiagonal(mol_ham_par);
    BOOST_CHECK(hubbard_diagonal.isApprox(fci_diagonal));

    Eigen::VectorXd hubbard_matvec = hubbard.matrixVectorProduct(mol_ham_par, hubbard_diagonal, hubbard_diagonal);
    Eigen::VectorXd fci_matvec = fci.matrixVectorProduct(mol_ham_par, fci_diagonal, fci_diagonal);
    BOOST_CHECK(hubbard_matvec.isApprox(fci_matvec));
}
