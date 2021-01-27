// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE "AdjacencyMatrix"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/ModelHamiltonian/AdjacencyMatrix.hpp"


/**
 *  Check if the adjacency matrix of a C3 graph is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(C3) {

    // Provide the reference value and check the result.
    auto A_ref = GQCP::SquareMatrix<size_t>::Zero(3);
    // clang-format off
    A_ref << 0, 1, 1,
             1, 0, 1,
             1, 1, 0;
    // clang-format on

    BOOST_CHECK(GQCP::AdjacencyMatrix::Cyclic(3).matrix().isApprox(A_ref));
}


/**
 *  Check if the adjacency matrix of a P4 graph is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(P4) {

    // Provide the reference value and check the result.
    auto A_ref = GQCP::SquareMatrix<size_t>::Zero(4);
    // clang-format off
    A_ref << 0, 1, 0, 0,
             1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0;
    // clang-format on

    BOOST_CHECK(GQCP::AdjacencyMatrix::Linear(4).matrix().isApprox(A_ref));
}


/**
 *  Check if the adjacency matrix of a P6 graph is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(P6) {

    // Provide the reference value and check the result.
    auto A_ref = GQCP::SquareMatrix<size_t>::Zero(6);
    // clang-format off
    A_ref << 0, 1, 0, 0, 0, 0,
             1, 0, 1, 0, 0, 0,
             0, 1, 0, 1, 0, 0,
             0, 0, 1, 0, 1, 0,
             0, 0, 0, 1, 0, 1,
             0, 0, 0, 0, 1, 0;
    // clang-format on

    BOOST_CHECK(GQCP::AdjacencyMatrix::Linear(6).matrix().isApprox(A_ref));
}


/**
 *  Check if the adjacency matrix of a P6 graph is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(C6) {

    // Provide the reference value and check the result.
    auto A_ref = GQCP::SquareMatrix<size_t>::Zero(6);
    // clang-format off
    A_ref << 0, 1, 0, 0, 0, 1,
             1, 0, 1, 0, 0, 0,
             0, 1, 0, 1, 0, 0,
             0, 0, 1, 0, 1, 0,
             0, 0, 0, 1, 0, 1,
             1, 0, 0, 0, 1, 0;
    // clang-format on

    BOOST_CHECK(GQCP::AdjacencyMatrix::Cyclic(6).matrix().isApprox(A_ref));
}
