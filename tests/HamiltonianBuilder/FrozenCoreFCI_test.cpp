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
#define BOOST_TEST_MODULE "FrozenCoreFCI"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "HamiltonianBuilder/FrozenCoreFCI.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/SelectedCI.hpp"




BOOST_AUTO_TEST_CASE ( FrozenCoreFCI_constructor ) {

    // Check if a correct constructor works
    GQCP::ProductFockSpace fock_space (15, 3, 3);
    GQCP::FrozenProductFockSpace frozen_fock_space (fock_space, 2);
    BOOST_CHECK_NO_THROW(GQCP::FrozenCoreFCI fci (frozen_fock_space));
}


BOOST_AUTO_TEST_CASE ( FrozenCoreFCI_public_methods ) {

    // Create random HamiltonianParameters to check compatibility
    size_t K = 5;
    auto random_hamiltonian_parameters = GQCP::HamiltonianParameters::Random(K);


    // Create a compatible Fock space
    GQCP::FrozenProductFockSpace fock_space (K, 3, 3, 1);
    GQCP::FrozenCoreFCI random_fci (fock_space);
    Eigen::VectorXd x = random_fci.calculateDiagonal(random_hamiltonian_parameters);
    BOOST_CHECK_NO_THROW(random_fci.constructHamiltonian(random_hamiltonian_parameters));
    BOOST_CHECK_NO_THROW(random_fci.matrixVectorProduct(random_hamiltonian_parameters, x, x));


    // Create an incompatible Fock space
    GQCP::FrozenProductFockSpace fock_space_invalid (K+1, 3, 3, 1);
    GQCP::FrozenCoreFCI random_fci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_fci_invalid.constructHamiltonian(random_hamiltonian_parameters), std::invalid_argument);
    BOOST_CHECK_THROW(random_fci_invalid.matrixVectorProduct(random_hamiltonian_parameters, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( SelectedCI_vs_FrozenCoreFCI ) {

    // Create H-chain HamiltonianParameters to test results
    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    auto random_hamiltonian_parameters = GQCP::HamiltonianParameters::Molecular(H5, "STO-3G");

    // Create compatible Fock spaces
    GQCP::FrozenProductFockSpace product_fock_space (K, 3, 3, 1);
    GQCP::SelectedFockSpace fock_space (product_fock_space);

    // The SelectedFockSpace includes the same configurations as the FrozenProductFockSpace
    // These builder instances should return the same results.
    GQCP::SelectedCI random_sci (fock_space);
    GQCP::FrozenCoreFCI random_frozen_core_fci (product_fock_space);

    Eigen::VectorXd sci_diagonal = random_sci.calculateDiagonal(random_hamiltonian_parameters);
    Eigen::VectorXd fci_diagonal = random_frozen_core_fci.calculateDiagonal(random_hamiltonian_parameters);

    Eigen::VectorXd sci_matvec = random_sci.matrixVectorProduct(random_hamiltonian_parameters, sci_diagonal, sci_diagonal);
    Eigen::VectorXd fci_matvec = random_frozen_core_fci.matrixVectorProduct(random_hamiltonian_parameters, fci_diagonal, fci_diagonal);

    Eigen::MatrixXd sci_ham = random_sci.constructHamiltonian(random_hamiltonian_parameters);
    Eigen::MatrixXd fci_ham = random_frozen_core_fci.constructHamiltonian(random_hamiltonian_parameters);

    BOOST_CHECK(sci_diagonal.isApprox(fci_diagonal));
    BOOST_CHECK(sci_matvec.isApprox(fci_matvec));
    BOOST_CHECK(sci_ham.isApprox(fci_ham));
}
