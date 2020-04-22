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

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"


/**
 *  Check the DOCI constructor.
 */
BOOST_AUTO_TEST_CASE(DOCI_constructor) {

    // Create a compatible ONV basis
    const GQCP::SeniorityZeroONVBasis onv_basis {15, 3};

    // Check if a correct constructor works
    BOOST_CHECK_NO_THROW(GQCP::DOCI {onv_basis});
}


/**
 *  Check the DOCI public methods.
 */
BOOST_AUTO_TEST_CASE(DOCI_public_methods) {

    // Create random HamiltonianParameters to check compatibility.
    const size_t K = 5;  // the number of spatial orbitals
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);

    // Create a compatible seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 3};
    const GQCP::DOCI doci_builder {onv_basis};

    // Test the public DOCI methods
    const auto x = doci_builder.calculateDiagonal(sq_hamiltonian);
    BOOST_CHECK_NO_THROW(doci_builder.constructHamiltonian(sq_hamiltonian));
    BOOST_CHECK_NO_THROW(doci_builder.matrixVectorProduct(sq_hamiltonian, x, x));

    // Create an incompatible seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis_incompatible {K + 1, 3};  // too many spatial orbitals
    const GQCP::DOCI doci_builder_incompatible {onv_basis_incompatible};
    BOOST_CHECK_THROW(doci_builder_incompatible.constructHamiltonian(sq_hamiltonian), std::invalid_argument);
    BOOST_CHECK_THROW(doci_builder_incompatible.matrixVectorProduct(sq_hamiltonian, x, x), std::invalid_argument);
}
