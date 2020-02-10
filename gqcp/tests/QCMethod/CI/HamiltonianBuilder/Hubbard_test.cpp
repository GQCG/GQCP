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

#include "QCMethod/CI/HamiltonianBuilder/Hubbard.hpp"

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"


BOOST_AUTO_TEST_CASE ( Hubbard_constructor ) {

    // Check if a correct constructor works
    GQCP::SpinResolvedONVBasis fock_space (15, 3, 3);
    BOOST_CHECK_NO_THROW(GQCP::Hubbard Hubbard (fock_space));
}


BOOST_AUTO_TEST_CASE ( Hubbard_public_methods ) {

    size_t K = 7;
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);


    // Create a compatible ONV basis
    GQCP::SpinResolvedONVBasis fock_space (K, 3, 3);
    GQCP::Hubbard random_Hubbard (fock_space);
    GQCP::VectorX<double> x = random_Hubbard.calculateDiagonal(sq_hamiltonian);
    BOOST_CHECK_NO_THROW(random_Hubbard.constructHamiltonian(sq_hamiltonian));
    BOOST_CHECK_NO_THROW(random_Hubbard.matrixVectorProduct(sq_hamiltonian, x, x));


    // Create an incompatible ONV basis
    GQCP::SpinResolvedONVBasis fock_space_invalid (K+1, 3, 3);

    // Create Hubbard module
    GQCP::Hubbard random_Hubbard_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_Hubbard_invalid.constructHamiltonian(sq_hamiltonian), std::invalid_argument);
    BOOST_CHECK_THROW(random_Hubbard_invalid.matrixVectorProduct(sq_hamiltonian, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI ) {

    // Create the Hubbard Hamiltonian
    size_t K = 4;
    size_t N = 2;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Hubbard(H);


    // Check if the Hamiltonian matrix is equal for FCI and Hubbard
    GQCP::SpinResolvedONVBasis fock_space (K, N, N);  // dim = 36
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    GQCP::SquareMatrix<double> hubbard_ham = hubbard.constructHamiltonian(sq_hamiltonian);
    GQCP::SquareMatrix<double> fci_ham = fci.constructHamiltonian(sq_hamiltonian);
    BOOST_CHECK(hubbard_ham.isApprox(fci_ham));
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_large ) {

    // Create the Hubbard Hamiltonian
    size_t K = 6;
    size_t N = 3;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Hubbard(H);


    // Check if the Hamiltonian matrix is equal for FCI and Hubbard
    GQCP::SpinResolvedONVBasis fock_space (K, N, N);  // dim = 400
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    GQCP::SquareMatrix<double> hubbard_ham = hubbard.constructHamiltonian(sq_hamiltonian);
    GQCP::SquareMatrix<double> fci_ham = fci.constructHamiltonian(sq_hamiltonian);
    BOOST_CHECK(hubbard_ham.isApprox(fci_ham));
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_matvec ) {

    // Create the Hubbard Hamiltonian
    size_t K = 4;
    size_t N = 2;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Hubbard(H);


    // Check if FCI and Hubbard have the same matvec for the Hubbard Hamiltonian
    GQCP::SpinResolvedONVBasis fock_space (K, N, N);  // dim = 36

    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    GQCP::VectorX<double> hubbard_diagonal = hubbard.calculateDiagonal(sq_hamiltonian);
    GQCP::VectorX<double> fci_diagonal = fci.calculateDiagonal(sq_hamiltonian);
    BOOST_CHECK(hubbard_diagonal.isApprox(fci_diagonal));

    GQCP::VectorX<double> hubbard_matvec = hubbard.matrixVectorProduct(sq_hamiltonian, hubbard_diagonal, hubbard_diagonal);
    GQCP::VectorX<double> fci_matvec = fci.matrixVectorProduct(sq_hamiltonian, fci_diagonal, fci_diagonal);
    BOOST_CHECK(hubbard_matvec.isApprox(fci_matvec));
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_large_matvec ) {

    // Create the Hubbard Hamiltonian
    size_t K = 6;
    size_t N = 3;
    auto H = GQCP::HoppingMatrix::Random(K);
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Hubbard(H);


    // Check if FCI and Hubbard have the same matvec for the Hubbard Hamiltonian
    GQCP::SpinResolvedONVBasis fock_space (K, N, N);  // dim = 400

    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    GQCP::VectorX<double> hubbard_diagonal = hubbard.calculateDiagonal(sq_hamiltonian);
    GQCP::VectorX<double> fci_diagonal = fci.calculateDiagonal(sq_hamiltonian);
    BOOST_CHECK(hubbard_diagonal.isApprox(fci_diagonal));

    GQCP::VectorX<double> hubbard_matvec = hubbard.matrixVectorProduct(sq_hamiltonian, hubbard_diagonal, hubbard_diagonal);
    GQCP::VectorX<double> fci_matvec = fci.matrixVectorProduct(sq_hamiltonian, fci_diagonal, fci_diagonal);
    BOOST_CHECK(hubbard_matvec.isApprox(fci_matvec));
}
