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

#define BOOST_TEST_MODULE "BlockMatrix_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/BlockMatrix.hpp"


BOOST_AUTO_TEST_CASE(constructor) {

    // Imagine the following 2x2 block is a part of of a 4x4 matrix:
    // x 1 2 x
    // x 3 4 x
    // x x x x
    // x x x x
    GQCP::MatrixX<size_t> block {2, 2};
    // clang-format off
    block << 1, 2,
             3, 4;
    // clang-format on

    // The block covers rows 0,2(not included) and columns 1,3(not included) of the full matrix
    GQCP::BlockMatrix<size_t> B1 {0, 2, 1, 3, block};
    BOOST_CHECK_THROW(GQCP::BlockMatrix<size_t>(0, 1, 1, 3, block), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE(operator_call) {

    // Imagine the following 2x2 block is a part of of a 4x4 matrix:
    // x 1 2 x
    // x 3 4 x
    // x x x x
    // x x x x
    GQCP::MatrixX<size_t> block {2, 2};
    // clang-format off
    block << 1, 2,
             3, 4;
    // clang-format on

    // The block covers rows 0,2(not included) and columns 1,3(not included) of the full matrix
    GQCP::BlockMatrix<size_t> B {0, 2, 1, 3, block};

    // The BlockMatrix B's operator() should behave like so:
    BOOST_CHECK(B(0, 1) == 1);
    BOOST_CHECK(B(0, 2) == 2);
    BOOST_CHECK(B(1, 1) == 3);
    BOOST_CHECK(B(1, 2) == 4);
}


/**
 *  Check the implementation of operator() for AP1roG-like variables, i.e. of the form G_i^a with 'i' an occupied orbital index and 'a' a virtual orbital index
 */
BOOST_AUTO_TEST_CASE(operator_call_AP1roG) {

    // Make an example for geminal coefficients for N_P=2 and K=5
    //      . .  1 2 3
    //      . .  4 5 6
    const size_t N_P = 2;  // the number of electron pairs
    const size_t K = 5;    // the number of spatial orbitals

    const size_t rows = N_P;      // the number of rows in the reduced geminal coefficient matrix
    const size_t cols = K - N_P;  // the number of columns in the reduced geminal coefficient matrix

    GQCP::VectorX<double> g {rows * cols};
    g << 1, 2, 3, 4, 5, 6;


    // Create the corresponding block matrix and check the implementation of operator()
    const GQCP::MatrixX<double> M = GQCP::MatrixX<double>::FromRowMajorVector(g, rows, cols);  // the actual geminal coefficients, reshaped into a matrix
    const GQCP::BlockMatrix<double> variables(0, N_P, N_P, K, M);                              // an encapsulating object that implements operator() intuitively. The block covers rows [0, N_P[ and columns [N_P, K[

    BOOST_CHECK_EQUAL(variables(0, 2), 1);
    BOOST_CHECK_EQUAL(variables(0, 3), 2);
    BOOST_CHECK_EQUAL(variables(0, 4), 3);
    BOOST_CHECK_EQUAL(variables(1, 2), 4);
    BOOST_CHECK_EQUAL(variables(1, 3), 5);
    BOOST_CHECK_EQUAL(variables(1, 4), 6);
}
