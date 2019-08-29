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
#define BOOST_TEST_MODULE "BlockMatrix_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/BlockMatrix.hpp"


BOOST_AUTO_TEST_CASE ( constructor ) {

    // Imagine the following 2x2 block is a part of of a 4x4 matrix:
    // x 1 2 x
    // x 3 4 x
    // x x x x
    // x x x x
    GQCP::MatrixX<size_t> block (2, 2);
    block << 1, 2,
             3, 4;

    // The block covers rows 0,2(not included) and columns 1,3(not included) of the full matrix
    GQCP::BlockMatrix<size_t> B1 (0,2, 1,3, block);
    BOOST_CHECK_THROW(GQCP::BlockMatrix<size_t>(0,1, 1,3, block), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( operator_call ) {

    // Imagine the following 2x2 block is a part of of a 4x4 matrix:
    // x 1 2 x
    // x 3 4 x
    // x x x x
    // x x x x
    GQCP::MatrixX<size_t> block (2, 2);
    block << 1, 2,
             3, 4;

    // The block covers rows 0,2(not included) and columns 1,3(not included) of the full matrix
    GQCP::BlockMatrix<size_t> B (0,2, 1,3, block);

    // The BlockMatrix B's operator() should behave like so:
    BOOST_CHECK(B(0,1) == 1);
    BOOST_CHECK(B(0,2) == 2);
    BOOST_CHECK(B(1,1) == 3);
    BOOST_CHECK(B(1,2) == 4);
}
