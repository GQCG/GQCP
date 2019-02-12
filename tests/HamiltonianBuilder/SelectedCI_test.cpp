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
#define BOOST_TEST_MODULE "SelectedCI"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "HamiltonianBuilder/SelectedCI.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianBuilder/DOCI.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"




BOOST_AUTO_TEST_CASE ( SelectedCI_constructor ) {

    // Check if a correct constructor works
    GQCP::ProductFockSpace product_fock_space (5, 3, 3);
    GQCP::SelectedFockSpace fock_space (product_fock_space);
    BOOST_CHECK_NO_THROW(GQCP::SelectedCI selected_ci (fock_space));
}


BOOST_AUTO_TEST_CASE ( SelectedCI_public_methods ) {

    // Create random HamiltonianParameters to check compatibility of the arguments
    size_t K = 5;
    auto random_hamiltonian_parameters = GQCP::HamiltonianParameters::Random(K);

    // Create a compatible Fock space
    GQCP::ProductFockSpace product_fock_space (K, 3, 3);
    GQCP::SelectedFockSpace fock_space (product_fock_space);
    GQCP::SelectedCI random_selected_ci (fock_space);
    Eigen::VectorXd x = random_selected_ci.calculateDiagonal(random_hamiltonian_parameters);
    BOOST_CHECK_NO_THROW(random_selected_ci.constructHamiltonian(random_hamiltonian_parameters));
    BOOST_CHECK_NO_THROW(random_selected_ci.matrixVectorProduct(random_hamiltonian_parameters, x, x));

    // Create an incompatible Fock space
    GQCP::ProductFockSpace product_fock_space_invalid (K+1, 3, 3);
    GQCP::SelectedFockSpace fock_space_invalid (product_fock_space_invalid);
    GQCP::SelectedCI random_selected_ci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_selected_ci_invalid.constructHamiltonian(random_hamiltonian_parameters), std::invalid_argument);
    BOOST_CHECK_THROW(random_selected_ci_invalid.matrixVectorProduct(random_hamiltonian_parameters, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( SelectedCI_vs_FCI ) {

    // Create H-chain HamiltonianParameters to test results from FCI and selected CI
    size_t K = 4;
    GQCP::Molecule H4 = GQCP::Molecule::HChain(K, 1.1);
    auto hamiltonian_parameters = GQCP::HamiltonianParameters::Molecular(H4, "STO-3G");

    // Create compatible Fock spaces
    GQCP::ProductFockSpace product_fock_space (K, 2, 2);
    GQCP::SelectedFockSpace fock_space (product_fock_space);

    // The SelectedFockSpace includes the same configurations as the ProductFockSpace, so the HamiltonianBuilders should return the same results
    GQCP::SelectedCI selected_ci (fock_space);
    GQCP::FCI fci (product_fock_space);

    Eigen::VectorXd selected_ci_diagonal = selected_ci.calculateDiagonal(hamiltonian_parameters);
    Eigen::VectorXd fci_diagonal = fci.calculateDiagonal(hamiltonian_parameters);

    Eigen::VectorXd selected_ci_matvec = selected_ci.matrixVectorProduct(hamiltonian_parameters, selected_ci_diagonal, selected_ci_diagonal);
    Eigen::VectorXd fci_matvec = fci.matrixVectorProduct(hamiltonian_parameters, fci_diagonal, fci_diagonal);

    Eigen::MatrixXd selected_ci_hamiltonian = selected_ci.constructHamiltonian(hamiltonian_parameters);
    Eigen::MatrixXd fci_hamiltonian = fci.constructHamiltonian(hamiltonian_parameters);

    BOOST_CHECK(selected_ci_diagonal.isApprox(fci_diagonal));
    BOOST_CHECK(selected_ci_matvec.isApprox(fci_matvec));
    BOOST_CHECK(selected_ci_hamiltonian.isApprox(fci_hamiltonian));
}

BOOST_AUTO_TEST_CASE ( SelectedCI_vs_DOCI ) {

    // Create H-chain HamiltonianParameters to test results from DOCI and selected CI
    size_t K = 4;
    GQCP::Molecule H4 = GQCP::Molecule::HChain(K, 1.1);
    auto hamiltonian_parameters = GQCP::HamiltonianParameters::Molecular(H4, "STO-3G");

    // Create compatible Fock spaces
    GQCP::FockSpace do_fock_space (K, 2);
    GQCP::SelectedFockSpace fock_space (do_fock_space);

    // The SelectedFockSpace includes the same configurations as the FockSpace
    // These builder instances should return the same results.
    GQCP::SelectedCI selected_ci (fock_space);
    GQCP::DOCI doci (do_fock_space);

    Eigen::VectorXd selected_ci_diagonal = selected_ci.calculateDiagonal(hamiltonian_parameters);
    Eigen::VectorXd doci_diagonal = doci.calculateDiagonal(hamiltonian_parameters);

    Eigen::VectorXd selected_ci_matvec = selected_ci.matrixVectorProduct(hamiltonian_parameters, selected_ci_diagonal, selected_ci_diagonal);
    Eigen::VectorXd doci_matvec = doci.matrixVectorProduct(hamiltonian_parameters, doci_diagonal, doci_diagonal);

    Eigen::MatrixXd selected_ci_hamiltonian = selected_ci.constructHamiltonian(hamiltonian_parameters);
    Eigen::MatrixXd doci_hamiltonian = doci.constructHamiltonian(hamiltonian_parameters);

    BOOST_CHECK(selected_ci_diagonal.isApprox(doci_diagonal));
    BOOST_CHECK(selected_ci_matvec.isApprox(doci_matvec));
    BOOST_CHECK(selected_ci_hamiltonian.isApprox(doci_hamiltonian));
}
