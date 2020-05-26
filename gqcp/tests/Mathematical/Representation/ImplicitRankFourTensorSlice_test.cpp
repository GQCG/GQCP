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

#define BOOST_TEST_MODULE "ImplicitRankFourTensorSlice_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/ImplicitRankFourTensorSlice.hpp"


/**
 *  Check if the call operator works as expected.
 */
BOOST_AUTO_TEST_CASE(operator_call) {

    // Create an implicit rank-four tensor and set some values through the call operator.
    auto B = GQCP::ImplicitRankFourTensorSlice<size_t>::ZeroFromBlockRanges(1, 3, 1, 3,
                                                                            1, 3, 1, 3);

    B(1, 1, 1, 1) = 1;
    B(1, 2, 2, 2) = 2;
    B(1, 2, 2, 1) = 3;
    B(1, 2, 1, 2) = 4;


    // Extract the updated dense slice representation and check if the values are filled in at the correct positions.
    const auto dense_slice_representation = B.asTensor();
    BOOST_CHECK(dense_slice_representation(0, 0, 0, 0) == 1);
    BOOST_CHECK(dense_slice_representation(0, 1, 1, 1) == 2);
    BOOST_CHECK(dense_slice_representation(0, 1, 1, 0) == 3);
    BOOST_CHECK(dense_slice_representation(0, 1, 0, 1) == 4);
}
