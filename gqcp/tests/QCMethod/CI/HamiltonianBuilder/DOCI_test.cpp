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
#define BOOST_TEST_MODULE "DOCI"

#include <boost/test/unit_test.hpp>

#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/RHF/PlainRHFSCFSolver.hpp"


BOOST_AUTO_TEST_CASE ( DOCI_constructor ) {
    // Create a compatible Fock space
    GQCP::FockSpace fock_space (15, 3);

    // Check if a correct constructor works
    BOOST_CHECK_NO_THROW(GQCP::DOCI doci (fock_space));
}


BOOST_AUTO_TEST_CASE ( DOCI_public_methods ) {

    // Create random HamiltonianParameters to check compatibility
    size_t K = 5;
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);

    // Create a compatible Fock space
    GQCP::FockSpace fock_space (K, 3);
    GQCP::DOCI random_doci (fock_space);

    // Test the public DOCI methods
    GQCP::VectorX<double> x = random_doci.calculateDiagonal(sq_hamiltonian);
    BOOST_CHECK_NO_THROW(random_doci.constructHamiltonian(sq_hamiltonian));
    BOOST_CHECK_NO_THROW(random_doci.matrixVectorProduct(sq_hamiltonian, x, x));

    // Create an incompatible Fock space
    GQCP::FockSpace fock_space_invalid (K+1, 3);
    GQCP::DOCI random_doci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_doci_invalid.constructHamiltonian(sq_hamiltonian), std::invalid_argument);
    BOOST_CHECK_THROW(random_doci_invalid.matrixVectorProduct(sq_hamiltonian, x, x), std::invalid_argument);
}
