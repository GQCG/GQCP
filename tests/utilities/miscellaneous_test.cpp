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
#define BOOST_TEST_MODULE "miscellaneous"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "utilities/miscellaneous.hpp"


BOOST_AUTO_TEST_CASE ( gray_code ) {

    BOOST_CHECK(GQCP::gray_code(0) == 0);    // "0000" (0)
    BOOST_CHECK(GQCP::gray_code(1) == 1);    // "0001" (1)
    BOOST_CHECK(GQCP::gray_code(2) == 3);    // "0011" (3)
    BOOST_CHECK(GQCP::gray_code(3) == 2);    // "0010" (2)
    BOOST_CHECK(GQCP::gray_code(4) == 6);    // "0110" (6)
    BOOST_CHECK(GQCP::gray_code(5) == 7);    // "0111" (7)
    BOOST_CHECK(GQCP::gray_code(6) == 5);    // "0101" (5)
    BOOST_CHECK(GQCP::gray_code(7) == 4);    // "0100" (4)
    BOOST_CHECK(GQCP::gray_code(8) == 12);   // "1100" (12)
    BOOST_CHECK(GQCP::gray_code(9) == 13);   // "1101" (13)
    BOOST_CHECK(GQCP::gray_code(10) == 15);  // "1111" (15)
    BOOST_CHECK(GQCP::gray_code(11) == 14);  // "1110" (14)
    BOOST_CHECK(GQCP::gray_code(12) == 10);  // "1010" (10)
    BOOST_CHECK(GQCP::gray_code(13) == 11);  // "1011" (11)
    BOOST_CHECK(GQCP::gray_code(14) == 9);   // "1001" (9)
    BOOST_CHECK(GQCP::gray_code(15) == 8);   // "1000" (8)
}


BOOST_AUTO_TEST_CASE ( vectorIndex ) {

    size_t cols = 11;
    size_t skipped = 2;
    BOOST_CHECK_EQUAL(GQCP::vectorIndex(0, 2, cols, skipped), 0);
    BOOST_CHECK_EQUAL(GQCP::vectorIndex(1, 2, cols, skipped), 9);

    cols = 5;
    // skipped = 0
    BOOST_CHECK_EQUAL(GQCP::vectorIndex(0, 2, cols), 2);
    BOOST_CHECK_EQUAL(GQCP::vectorIndex(1, 2, cols), 7);
}


BOOST_AUTO_TEST_CASE ( matrixIndex ) {

    size_t cols = 11;
    size_t skipped = 2;

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
