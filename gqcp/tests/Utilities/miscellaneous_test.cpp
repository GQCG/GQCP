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

#define BOOST_TEST_MODULE "miscellaneous"

#include <boost/test/unit_test.hpp>

#include "Utilities/miscellaneous.hpp"


/**
 *  Check if the Gray code is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(grayCodeOf) {

    BOOST_CHECK(GQCP::grayCodeOf(0) == 0);    // "0000" (0)
    BOOST_CHECK(GQCP::grayCodeOf(1) == 1);    // "0001" (1)
    BOOST_CHECK(GQCP::grayCodeOf(2) == 3);    // "0011" (3)
    BOOST_CHECK(GQCP::grayCodeOf(3) == 2);    // "0010" (2)
    BOOST_CHECK(GQCP::grayCodeOf(4) == 6);    // "0110" (6)
    BOOST_CHECK(GQCP::grayCodeOf(5) == 7);    // "0111" (7)
    BOOST_CHECK(GQCP::grayCodeOf(6) == 5);    // "0101" (5)
    BOOST_CHECK(GQCP::grayCodeOf(7) == 4);    // "0100" (4)
    BOOST_CHECK(GQCP::grayCodeOf(8) == 12);   // "1100" (12)
    BOOST_CHECK(GQCP::grayCodeOf(9) == 13);   // "1101" (13)
    BOOST_CHECK(GQCP::grayCodeOf(10) == 15);  // "1111" (15)
    BOOST_CHECK(GQCP::grayCodeOf(11) == 14);  // "1110" (14)
    BOOST_CHECK(GQCP::grayCodeOf(12) == 10);  // "1010" (10)
    BOOST_CHECK(GQCP::grayCodeOf(13) == 11);  // "1011" (11)
    BOOST_CHECK(GQCP::grayCodeOf(14) == 9);   // "1001" (9)
    BOOST_CHECK(GQCP::grayCodeOf(15) == 8);   // "1000" (8)
}


/**
 *  Check if `vectorIndex` behaves correctly.
 */
BOOST_AUTO_TEST_CASE(vectorIndex) {

    size_t cols = 11;
    const size_t skipped = 2;
    BOOST_CHECK_EQUAL(GQCP::vectorIndex(0, 2, cols, skipped), 0);
    BOOST_CHECK_EQUAL(GQCP::vectorIndex(1, 2, cols, skipped), 9);

    cols = 5;
    // skipped = 0
    BOOST_CHECK_EQUAL(GQCP::vectorIndex(0, 2, cols), 2);
    BOOST_CHECK_EQUAL(GQCP::vectorIndex(1, 2, cols), 7);
}


/**
 *  Check if `matrixIndex` behaves correctly.
 */
BOOST_AUTO_TEST_CASE(matrixIndex) {

    size_t cols = 11;
    const size_t skipped = 2;

    BOOST_CHECK_EQUAL(GQCP::matrixIndexMajor(0, cols, skipped), 0);
    BOOST_CHECK_EQUAL(GQCP::matrixIndexMajor(9, cols, skipped), 1);

    BOOST_CHECK_EQUAL(GQCP::matrixIndexMinor(0, cols, skipped), 2);
    BOOST_CHECK_EQUAL(GQCP::matrixIndexMinor(4, cols, skipped), 6);


    cols = 5;
    // skipped = 0
    BOOST_CHECK_EQUAL(GQCP::matrixIndexMajor(2, cols), 0);
    BOOST_CHECK_EQUAL(GQCP::matrixIndexMajor(6, cols), 1);

    BOOST_CHECK_EQUAL(GQCP::matrixIndexMinor(2, cols), 2);
    BOOST_CHECK_EQUAL(GQCP::matrixIndexMinor(7, cols), 2);
}


/**
 *  Check if generatePartitionsOf works as expected.
 */
BOOST_AUTO_TEST_CASE(generatePartitionsOf) {

    // 2-way partition '2'.
    const std::vector<std::vector<size_t>> ref_partitions1 {{2, 0}, {1, 1}, {0, 2}};
    BOOST_CHECK(GQCP::generatePartitionsOf(2, 2) == ref_partitions1);


    // 3-way partition '2'.
    const std::vector<std::vector<size_t>> ref_partitions2 {{2, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 2, 0}, {0, 1, 1}, {0, 0, 2}};
    BOOST_CHECK(GQCP::generatePartitionsOf(2, 3) == ref_partitions2);


    // 3-way partition '3'.
    const std::vector<std::vector<size_t>> ref_partitions3 {{3, 0, 0}, {2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 1, 1}, {1, 0, 2}, {0, 3, 0}, {0, 2, 1}, {0, 1, 2}, {0, 0, 3}};
    BOOST_CHECK(GQCP::generatePartitionsOf(3, 3) == ref_partitions3);


    // 4-way partition '2'.
    const std::vector<std::vector<size_t>> ref_partitions4 {{2, 0, 0, 0}, {1, 1, 0, 0}, {1, 0, 1, 0}, {1, 0, 0, 1}, {0, 2, 0, 0}, {0, 1, 1, 0}, {0, 1, 0, 1}, {0, 0, 2, 0}, {0, 0, 1, 1}, {0, 0, 0, 2}};
    BOOST_CHECK(GQCP::generatePartitionsOf(2, 4) == ref_partitions4);
}


/**
 *  Check if generateUniquePartitionsOf works as expected.
 */
BOOST_AUTO_TEST_CASE(generateUniquePartitionsOf) {

    // 2-way partition '4'.
    const std::vector<std::vector<size_t>> ref_partitions1 {{4, 0}, {3, 1}, {2, 2}};
    BOOST_CHECK(GQCP::generateUniquePartitionsOf(4, 2) == ref_partitions1);


    // 3-way partition '2'.
    const std::vector<std::vector<size_t>> ref_partitions2 {{2, 0, 0}, {1, 1, 0}};
    BOOST_CHECK(GQCP::generateUniquePartitionsOf(2, 3) == ref_partitions2);


    // 3-way partition '3'.
    const std::vector<std::vector<size_t>> ref_partitions3 {{3, 0, 0}, {2, 1, 0}, {1, 1, 1}};
    BOOST_CHECK(GQCP::generateUniquePartitionsOf(3, 3) == ref_partitions3);


    // 5-way partition '5'.
    const std::vector<std::vector<size_t>> ref_partitions4 {{5, 0, 0, 0, 0},
                                                            {4, 1, 0, 0, 0},
                                                            {3, 2, 0, 0, 0},
                                                            {3, 1, 1, 0, 0},
                                                            {2, 2, 1, 0, 0},
                                                            {2, 1, 1, 1, 0},
                                                            {1, 1, 1, 1, 1}};
    BOOST_CHECK(GQCP::generateUniquePartitionsOf(5, 5) == ref_partitions4);
}


/**
 *  Check the implementation of `triangularRoot` and `strictTriangularRoot`.
 */
BOOST_AUTO_TEST_CASE(triangularRoot_strictTriangularRoot) {

    BOOST_CHECK(GQCP::triangularRootOf(6) == 3);
    BOOST_CHECK(GQCP::strictTriangularRootOf(3) == 3);

    BOOST_CHECK(GQCP::triangularRootOf(10) == 4);
    BOOST_CHECK(GQCP::strictTriangularRootOf(6) == 4);
}


/**
 *  Check if `validateAndOpen` handles errors correctly.
 */
BOOST_AUTO_TEST_CASE(validateAndOpen) {

    // Make sure we get an error when a nonsense path is given (i.e. no extension).
    BOOST_REQUIRE_THROW(GQCP::validateAndOpen("this is a nonsense data path", "data"), std::invalid_argument);

    // Make sure we get an error when a path with a wrong extension is given.
    BOOST_REQUIRE_THROW(GQCP::validateAndOpen("data/small_vector.data", "xyz"), std::invalid_argument);

    // Make sure we don't get an error when a correct path is given.
    BOOST_REQUIRE_NO_THROW(GQCP::validateAndOpen("data/h2o.xyz", "xyz"));
}


/**
 *  Check if `findElementIndex` is correct.
 */
BOOST_AUTO_TEST_CASE(findElementIndex) {

    const std::vector<int> vector {1, 2, 3};

    BOOST_REQUIRE_THROW(GQCP::findElementIndex(vector, 0), std::out_of_range);  // 0 is not in the vector.
    BOOST_CHECK_EQUAL(GQCP::findElementIndex(vector, 2), 1);
}
