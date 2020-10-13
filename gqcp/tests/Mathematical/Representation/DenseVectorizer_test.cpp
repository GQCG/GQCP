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

#define BOOST_TEST_MODULE "DenseVectorizer_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/DenseVectorizer.hpp"


/**
 *  Check if the offset is calculated correctly when rowMajor ordering is used.
 */
BOOST_AUTO_TEST_CASE(offset_row) {

    // To understand the offset method, let's take a look at an example using a 3x3 matrix

    // [[1, 2, 3],
    //  [4, 5, 6],
    //  [7, 8, 9]]

    // Using row-major order, this 2D array can be seen as:
    // [1, 2, 3, 4, 5, 6, 7, 8, 9]

    // The strides for this matrix are [3, 1]
    // To move to the next row you move three positions in the 1D array, to move to the next column you move one position

    // Create an array containing the dimension of each axis
    // Create a vectorizer with axis dimensions 2, thus representing a 2D array and define its axis dimensions with the previously created array
    const auto dimensions = std::array<size_t, 2> {3, 3};
    const auto test_vectorizer = GQCP::DenseVectorizer<2>(dimensions, GQCP::Ordering::RowMajor);

    // The offset method only works for indices smaller than the axis dimension, since you can't move to an index greater than the axis dimensions
    // Check whther an x-index being too large causes a throw
    const std::array<size_t, 2> indices_wrong_1 = {4, 2};
    BOOST_CHECK_THROW(test_vectorizer.offset(indices_wrong_1), std::invalid_argument);

    // Check whether a y-index being too large causes a throw
    const std::array<size_t, 2> indices_wrong_2 = {2, 4};
    BOOST_CHECK_THROW(test_vectorizer.offset(indices_wrong_2), std::invalid_argument);

    // Move one column and check the offset
    const int ref_1 = 1;
    const std::array<size_t, 2> indices_1 = {0, 1};
    BOOST_CHECK(test_vectorizer.offset(indices_1) == ref_1);

    // Move one row and check the offset
    const int ref_2 = 3;
    const std::array<size_t, 2> indices_2 = {1, 0};
    BOOST_CHECK(test_vectorizer.offset(indices_2) == ref_2);
}


/**
 *  Check if the offset is calculated correctly when columnMajor ordering is used.
 */
BOOST_AUTO_TEST_CASE(offset_column) {

    // To understand the offset method, let's take a look at an example using a 3x3 matrix

    // [[1, 2, 3],
    //  [4, 5, 6],
    //  [7, 8, 9]]

    // Using column-major order, this 2D array can be seen as:
    // [1, 4, 7, 2, 5, 8, 3, 6, 9]

    // The strides for this matrix are [1, 3]
    // To move to the next row you move one position in the 1D array, to move to the next column you move three positions

    // Create an array containing the dimension of each axis
    // Create a vectorizer with axis dimensions 2, thus representing a 2D array and define its axis dimensions with the previously created array
    const auto dimensions = std::array<size_t, 2> {3, 3};
    const auto test_vectorizer = GQCP::DenseVectorizer<2>(dimensions);

    // Move one column and check the offset
    const int ref_1 = 3;
    const std::array<size_t, 2> indices_1 = {0, 1};
    BOOST_CHECK(test_vectorizer.offset(indices_1) == ref_1);

    // Move one row and check the offset
    const int ref_2 = 1;
    const std::array<size_t, 2> indices_2 = {1, 0};
    BOOST_CHECK(test_vectorizer.offset(indices_2) == ref_2);
}

/**
 *  Check if the offset is calculated correctly when columnMajor ordering is used.
 */
BOOST_AUTO_TEST_CASE(offset_multi_dimensional) {

    // The offset method can also be used on multi dimensional arrays, let's take a look at an example using a 3D array.

    //  [[[ 0,  1,  2,  3],
    //    [ 4,  5,  6,  7],
    //    [ 8,  9, 10, 11]],
    //   [[12, 13, 14, 15],
    //    [16, 17, 18, 19],
    //    [20, 21, 22, 23]]]

    // Using column-major order, this 3D array can be seen as:
    // [0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11, 12, 16, 20, 13, 17, 21, 14, 18, 22, 15, 19, 23]

    // The strides for this matrix are [12, 1, 4]
    // To move to the next row you move one position in the 1D-array, to move to the next column you move three positions and to move to the next inner array you move twelve positions

    // Create an array containing the dimension of each axis
    // Create a vectorizer with axis dimensions 3, thus representing a 3D array and define its axis dimensions with the previously created array
    const auto dimensions = std::array<size_t, 3> {3, 4, 2};
    const auto test_vectorizer = GQCP::DenseVectorizer<3>(dimensions);

    // Move one inner array and check the offset
    const int ref_1 = 12;
    const std::array<size_t, 3> indices_1 = {0, 0, 1};
    BOOST_CHECK(test_vectorizer.offset(indices_1) == ref_1);

    // Move one column and two rows and check the offset
    const int ref_2 = 5;
    const std::array<size_t, 3> indices_2 = {2, 1, 0};
    BOOST_CHECK(test_vectorizer.offset(indices_2) == ref_2);

    // Move one row and three columns check the offset
    const int ref_3 = 10;
    const std::array<size_t, 3> indices_3 = {1, 3, 0};
    BOOST_CHECK(test_vectorizer.offset(indices_3) == ref_3);
}
