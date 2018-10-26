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
#define BOOST_TEST_MODULE "Hubbard"


#include "HamiltonianBuilder/Hubbard.hpp"

#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( Hubbard_constructor ) {
    // Create a compatible Fock space
    GQCP::ProductFockSpace fock_space (15, 3, 3);

    // Check if a correct constructor works
    BOOST_CHECK_NO_THROW(GQCP::Hubbard Hubbard (fock_space));
}


BOOST_AUTO_TEST_CASE ( Hubbard_public_methods ) {
    // Create an AOBasis
    GQCP::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(water, "STO-3G");


    // Create random HamiltonianParameters from One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis->get_number_of_basis_functions();
    GQCP::HamiltonianParameters random_hamiltonian_parameters = GQCP::constructRandomHamiltonianParameters(K);

    // Create a compatible Fock space
    GQCP::ProductFockSpace fock_space (K, 3, 3);

    // Create Hubbard module
    GQCP::Hubbard random_Hubbard (fock_space);

    // Test the public Hubbard methods
    Eigen::VectorXd x = random_Hubbard.calculateDiagonal(random_hamiltonian_parameters);
    BOOST_CHECK_NO_THROW(random_Hubbard.constructHamiltonian(random_hamiltonian_parameters));
    BOOST_CHECK_NO_THROW(random_Hubbard.matrixVectorProduct(random_hamiltonian_parameters, x, x));

    // Create an incompatible Fock space
    GQCP::ProductFockSpace fock_space_i (K+1, 3, 3);

    // Create Hubbard module
    GQCP::Hubbard random_Hubbard_i (fock_space_i);
    BOOST_CHECK_THROW(random_Hubbard_i.constructHamiltonian(random_hamiltonian_parameters), std::invalid_argument);
    BOOST_CHECK_THROW(random_Hubbard_i.matrixVectorProduct(random_hamiltonian_parameters, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI ) {

    // Check if FCI and Hubbard produce the Hamiltonian matrix for Hubbard Hamiltonian parameters

    // Create the Hamiltonian parameters for the triagonal of a Hubbard lattice.
    Eigen::VectorXd triagonal_test = Eigen::VectorXd::Random(10);

    size_t N = 2;
    auto mol_ham_par = GQCP::constructHubbardParameters(triagonal_test);
    auto K = mol_ham_par.get_K();


    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 36

    // Create the Hubbard module
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    Eigen::MatrixXd hubbard_ham = hubbard.constructHamiltonian(mol_ham_par);
    Eigen::MatrixXd fci_ham = fci.constructHamiltonian(mol_ham_par);

    BOOST_CHECK(hubbard_ham.isApprox(fci_ham));
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_large ) {

    // Check if FCI and Hubbard produce the Hamiltonian matrix for Hubbard Hamiltonian parameters

    // Create the Hamiltonian parameters for the triagonal of a Hubbard lattice.
    Eigen::VectorXd triagonal_test = Eigen::VectorXd::Random(21);

    size_t N = 3;
    auto mol_ham_par = GQCP::constructHubbardParameters(triagonal_test);
    auto K = mol_ham_par.get_K();


    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 400

    // Create the Hubbard module
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    Eigen::MatrixXd hubbard_ham = hubbard.constructHamiltonian(mol_ham_par);
    Eigen::MatrixXd fci_ham = fci.constructHamiltonian(mol_ham_par);

    BOOST_CHECK(hubbard_ham.isApprox(fci_ham));
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_matvec ) {

    // Check if FCI and Hubbard have the same matvec for Hubbard Hamiltonian parameters

    // Create the Hamiltonian parameters for the triagonal of a Hubbard lattice.
    Eigen::VectorXd triagonal_test = Eigen::VectorXd::Random(10);

    size_t N = 2;
    auto mol_ham_par = GQCP::constructHubbardParameters(triagonal_test);
    auto K = mol_ham_par.get_K();


    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 36

    // Create the Hubbard module
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    Eigen::VectorXd hubbard_diagonal = hubbard.calculateDiagonal(mol_ham_par);
    Eigen::VectorXd fci_diagonal = fci.calculateDiagonal(mol_ham_par);

    Eigen::VectorXd hubbard_matvec = hubbard.matrixVectorProduct(mol_ham_par, hubbard_diagonal, hubbard_diagonal);
    Eigen::VectorXd fci_matvec = fci.matrixVectorProduct(mol_ham_par, fci_diagonal, fci_diagonal);

    BOOST_CHECK(hubbard_matvec.isApprox(fci_matvec));
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_large_matvec ) {

    // Check if FCI and Hubbard have the same matvec for Hubbard Hamiltonian parameters

    // Create the Hamiltonian parameters for the triagonal of a Hubbard lattice.
    Eigen::VectorXd triagonal_test = Eigen::VectorXd::Random(21);

    size_t N = 3;
    auto mol_ham_par = GQCP::constructHubbardParameters(triagonal_test);
    auto K = mol_ham_par.get_K();


    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 400

    // Create the Hubbard module
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    Eigen::VectorXd hubbard_diagonal = hubbard.calculateDiagonal(mol_ham_par);
    Eigen::VectorXd fci_diagonal = fci.calculateDiagonal(mol_ham_par);

    Eigen::VectorXd hubbard_matvec = hubbard.matrixVectorProduct(mol_ham_par, hubbard_diagonal, hubbard_diagonal);
    Eigen::VectorXd fci_matvec = fci.matrixVectorProduct(mol_ham_par, fci_diagonal, fci_diagonal);

    BOOST_CHECK(hubbard_matvec.isApprox(fci_matvec));
}