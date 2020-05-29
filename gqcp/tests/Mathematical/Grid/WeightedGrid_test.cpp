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

#define BOOST_TEST_MODULE "WeightedGrid_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Grid/WeightedGrid.hpp"


/**
 *  Test if the constructor of WeightedGrid throws as expected.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    // Prepare some valid and invalid combinations of weights and points.
    GQCP::ArrayX<double> weights1 {3};
    weights1 << 1.0, 2.0, 3.0;

    std::vector<GQCP::Vector<double, 3>> points1;
    points1.emplace_back(1.0, 0.0, 0.0);
    points1.emplace_back(0.0, 1.0, 0.0);
    points1.emplace_back(0.0, 0.0, 1.0);

    GQCP::ArrayX<double> weights2 {4};
    weights2 << 1.0, 2.0, 3.0, 4.0;

    std::vector<GQCP::Vector<double, 3>> points2;
    points2.emplace_back(0.0, 0.0, 0.0);
    points2.emplace_back(1.0, 0.0, 0.0);
    points2.emplace_back(0.0, 1.0, 0.0);
    points2.emplace_back(0.0, 0.0, 1.0);


    // Check if the constructor throws as expected.
    BOOST_CHECK_NO_THROW(GQCP::WeightedGrid(points1, weights1));
    BOOST_CHECK_NO_THROW(GQCP::WeightedGrid(points2, weights2));

    BOOST_CHECK_THROW(GQCP::WeightedGrid(points1, weights2), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::WeightedGrid(points2, weights1), std::invalid_argument);
}
