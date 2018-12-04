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
#define BOOST_TEST_MODULE "FCI"


#include "HamiltonianBuilder/FCI.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( FCI_constructor ) {

    // Check if a correct constructor works
    GQCP::ProductFockSpace fock_space (15, 3, 3);
    BOOST_CHECK_NO_THROW(GQCP::FCI fci (fock_space));
}


BOOST_AUTO_TEST_CASE ( FCI_public_methods ) {

    // Create random HamiltonianParameters to check compatibility
    size_t K = 5;
    auto random_hamiltonian_parameters = GQCP::HamiltonianParameters::Random(K);


    // Create a compatible Fock space
    GQCP::ProductFockSpace fock_space (K, 3, 3);
    GQCP::FCI random_fci (fock_space);
    Eigen::VectorXd x = random_fci.calculateDiagonal(random_hamiltonian_parameters);
    BOOST_CHECK_NO_THROW(random_fci.constructHamiltonian(random_hamiltonian_parameters));
    BOOST_CHECK_NO_THROW(random_fci.matrixVectorProduct(random_hamiltonian_parameters, x, x));


    // Create an incompatible Fock space
    GQCP::ProductFockSpace fock_space_invalid (K+1, 3, 3);
    GQCP::FCI random_fci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_fci_invalid.constructHamiltonian(random_hamiltonian_parameters), std::invalid_argument);
    BOOST_CHECK_THROW(random_fci_invalid.matrixVectorProduct(random_hamiltonian_parameters, x, x), std::invalid_argument);
}
