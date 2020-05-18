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

#define BOOST_TEST_MODULE "ImplicitMatrixSlice_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/ImplicitMatrixSlice.hpp"


/**
 *  Check if the basic constructor works and throws as expected.
 */
BOOST_AUTO_TEST_CASE(basic_constructor) {

    // Imagine the following 2x2 slice is part of an implicit 4x4 matrix:
    // x 1 x 2
    // x 3 x 4
    // x x x x
    // x x x x

    // The dense matrix representation of the slice is then given by:
    GQCP::MatrixX<size_t> slice {2, 2};
    // clang-format off
    slice << 1, 2,
             3, 4;
    // clang-format on


    // The following unordered_maps map the implicit indices to the dense indices.
    const std::unordered_map<size_t, size_t> rows_map {
        {0, 0},  // the 0'th implicit row maps to the 0'th dense row
        {1, 1}   // the 1'th implicit row maps to the 1'th dense row
    };

    const std::unordered_map<size_t, size_t> cols_map {
        {1, 0},  // the 1'th implicit column maps to the 0'th dense column
        {3, 1}   // the 1'th implicit column maps to the 1'th dense column
    };


    BOOST_CHECK_NO_THROW(GQCP::ImplicitMatrixSlice<size_t>(rows_map, cols_map, slice));


    // Proceed to 'invalidate' the rows and columns mappings and check if the constructor throws as expected.
    auto rows_map_copy = rows_map;
    rows_map_copy[5] = 5;

    auto cols_map_copy = cols_map;
    cols_map_copy[5] = 5;


    BOOST_CHECK_THROW(GQCP::ImplicitMatrixSlice<size_t>(rows_map_copy, cols_map, slice), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::ImplicitMatrixSlice<size_t>(rows_map, cols_map_copy, slice), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::ImplicitMatrixSlice<size_t>(rows_map_copy, cols_map_copy, slice), std::invalid_argument);
}


/**
 *  Check if a constructor works as expected.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    // Imagine the following 2x2 block is a part of an implicit 4x4 matrix:
    // x 1 2 x
    // x 3 4 x
    // x x x x
    // x x x x
    GQCP::MatrixX<size_t> slice {2, 2};
    // clang-format off
    slice << 1, 2,
             3, 4;
    // clang-format on

    // The block covers rows 0,2(not included) and columns 1,3(not included) of the implicit matrix.
    GQCP::ImplicitMatrixSlice<size_t> B1 {0, 2, 1, 3, slice};                                        // this constructor should work
    BOOST_CHECK_THROW(GQCP::ImplicitMatrixSlice<size_t>(0, 1, 1, 3, slice), std::invalid_argument);  // this constructor should fail
}


/**
 *  Check if the call operator works as expected.
 */
BOOST_AUTO_TEST_CASE(operator_call) {

    // Imagine the following 2x2 block is a part of an implicit 4x4 matrix:
    // x 1 2 x
    // x 3 4 x
    // x x x x
    // x x x x
    GQCP::MatrixX<size_t> slice {2, 2};
    // clang-format off
    slice << 1, 2,
             3, 4;
    // clang-format on

    // The block covers rows 0,2(not included) and columns 1,3(not included) of the implicit matrix.
    GQCP::ImplicitMatrixSlice<size_t> B {0, 2, 1, 3, slice};

    // Check if ImplicitMatrixSlice B's operator() behaves as expected.
    BOOST_CHECK(B(0, 1) == 1);
    BOOST_CHECK(B(0, 2) == 2);
    BOOST_CHECK(B(1, 1) == 3);
    BOOST_CHECK(B(1, 2) == 4);
}


/**
 *  Check the implementation of operator() for occupied-virtual-like objects, like the AP1roG geminal coefficients.
 */
BOOST_AUTO_TEST_CASE(operator_call_AP1roG) {

    // Make an example for geminal coefficients for N_P=2 and K=5.
    //      . .  1 2 3
    //      . .  4 5 6
    const size_t N_P = 2;  // the number of electron pairs
    const size_t K = 5;    // the number of spatial orbitals

    const size_t rows = N_P;      // the number of rows in the reduced geminal coefficient matrix
    const size_t cols = K - N_P;  // the number of columns in the reduced geminal coefficient matrix

    GQCP::VectorX<double> g {rows * cols};
    g << 1, 2, 3, 4, 5, 6;


    // Create the corresponding slice and check the implementation of operator().
    const GQCP::MatrixX<double> M = GQCP::MatrixX<double>::FromRowMajorVector(g, rows, cols);  // the actual geminal coefficients, reshaped into a matrix
    const GQCP::ImplicitMatrixSlice<double> variables {0, N_P, N_P, K, M};                     // an encapsulating object that implements operator() intuitively. The block covers rows [0, N_P[ and columns [N_P, K[

    BOOST_CHECK_EQUAL(variables(0, 2), 1);
    BOOST_CHECK_EQUAL(variables(0, 3), 2);
    BOOST_CHECK_EQUAL(variables(0, 4), 3);
    BOOST_CHECK_EQUAL(variables(1, 2), 4);
    BOOST_CHECK_EQUAL(variables(1, 3), 5);
    BOOST_CHECK_EQUAL(variables(1, 4), 6);
}
