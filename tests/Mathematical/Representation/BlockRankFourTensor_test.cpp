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
#define BOOST_TEST_MODULE "BlockRankFourTensor_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/BlockRankFourTensor.hpp"


BOOST_AUTO_TEST_CASE ( operator_call ) {

    // A zero block tensor of dimension 2
    GQCP::BlockRankFourTensor<size_t> B (1,3, 1,3, 1,3, 1,3);


    // The BlockRankFourTensor B's operator() should behave like so:
    // Set some values
    B(1,1,1,1) = 1;
    B(1,2,2,2) = 2;
    B(1,2,2,1) = 3;
    B(1,2,1,2) = 4;

    // Extract the updated block and check if the values are filled in at the correct positions
    const auto block = B.asTensor();
    BOOST_CHECK(block(0,0,0,0) == 1);
    BOOST_CHECK(block(0,1,1,1) == 2);
    BOOST_CHECK(block(0,1,1,0) == 3);
    BOOST_CHECK(block(0,1,0,1) == 4);
}
