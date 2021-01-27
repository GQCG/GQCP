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

#define BOOST_TEST_MODULE "HoppingMatrix"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/ModelHamiltonian/HoppingMatrix.hpp"


/**
 *  Test if the HoppingMatrix constructor throws (and doesn't) as expected.
 */
BOOST_AUTO_TEST_CASE(constructor_throws) {

    GQCP::MatrixX<double> H1(3, 4);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix<double> H(H1), std::invalid_argument);  // Not square.

    GQCP::MatrixX<double> H2 = GQCP::MatrixX<double>::Random(3, 3);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix<double> H(H2), std::invalid_argument);  // Not symmetric.

    GQCP::MatrixX<double> H3 = H2 + H2.transpose();
    BOOST_CHECK_NO_THROW(GQCP::HoppingMatrix<double> H(H3));  // Square and symmetric: OK.
}


/**
 *  Check if the conversion from an adjacency matrix (and the parameters t and U) to a hopping matrix works as expected.
 */
BOOST_AUTO_TEST_CASE(triangle_adjacency_matrix) {

    const double t = 1.0;
    const double U = 2.0;

    GQCP::SquareMatrix<double> H_ref = GQCP::SquareMatrix<double>::Zero(3);
    // clang-format off
    H_ref << U, -t, -t,
            -t,  U, -t,
            -t, -t,  U;
    // clang-format on


    // Construct an adjacency matrix, convert it to a hopping matrix and check the result.
    const auto A = GQCP::AdjacencyMatrix::Cyclic(3);

    const GQCP::HoppingMatrix<double> H {A, t, U};
    BOOST_CHECK(H_ref.isApprox(H.matrix()));
}


/**
 *  Check if the conversion from a comma-separated line (CSLine) works as expected.
 */
BOOST_AUTO_TEST_CASE(FromCSLine) {

    GQCP::SquareMatrix<double> H_ref = GQCP::SquareMatrix<double>::Zero(3);
    // clang-format off
    H_ref << 1, 2, 4,
             2, 3, 5,
             4, 5, 6;
    // clang-format on

    const auto H = GQCP::HoppingMatrix<double>::FromCSLine("1,2,4,3,5,6");

    BOOST_CHECK(H_ref.isApprox(H.matrix()));
}
