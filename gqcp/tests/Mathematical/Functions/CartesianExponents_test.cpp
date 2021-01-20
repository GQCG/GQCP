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

#define BOOST_TEST_MODULE "CartesianExponents"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Functions/CartesianExponents.hpp"


/**
 *  Check if the angular momentum of CartesianExponents is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(angular_momentum) {

    const GQCP::CartesianExponents exponents {2, 3, 4};

    BOOST_CHECK_EQUAL(exponents.angularMomentum(), 9);
}


/**
 *  Check if `operator<` behaves as expected.
 */
BOOST_AUTO_TEST_CASE(operator_smaller) {

    const GQCP::CartesianExponents exponents1 {1, 0, 0};
    const GQCP::CartesianExponents exponents2 {2, 0, 0};
    const GQCP::CartesianExponents exponents3 {1, 1, 0};
    const GQCP::CartesianExponents exponents4 {0, 1, 1};

    BOOST_CHECK(exponents1 < exponents2);
    BOOST_CHECK(exponents1 < exponents3);
    BOOST_CHECK(exponents1 < exponents4);

    BOOST_CHECK(exponents2 < exponents3);
    BOOST_CHECK(exponents2 < exponents4);
    BOOST_CHECK(exponents3 < exponents4);
}


/**
 *  Check if `operator==` behaves as expected.
 */
BOOST_AUTO_TEST_CASE(operator_equals) {

    GQCP::CartesianExponents exponents1 {1, 0, 0};
    GQCP::CartesianExponents exponents2 {1, 0, 0};
    GQCP::CartesianExponents exponents3 {2, 0, 0};

    BOOST_CHECK(exponents1 == exponents2);
    BOOST_CHECK(!(exponents1 == exponents3));
}


/**
 *  Check if `allPermutations` generates the correct permutations of the Cartesian exponents.
 */
BOOST_AUTO_TEST_CASE(allPermutations) {

    // Example 1.
    const GQCP::CartesianExponents exponents1 {2, 3, 4};
    const auto all_permutations1 = exponents1.allPermutations();
    const std::vector<GQCP::CartesianExponents> ref_all_permutations1 {
        GQCP::CartesianExponents(4, 3, 2), GQCP::CartesianExponents(4, 2, 3),
        GQCP::CartesianExponents(3, 4, 2), GQCP::CartesianExponents(3, 2, 4),
        GQCP::CartesianExponents(2, 4, 3), GQCP::CartesianExponents(2, 3, 4)};
    BOOST_CHECK(all_permutations1 == ref_all_permutations1);


    // Example 2.
    const GQCP::CartesianExponents exponents2 {1, 0, 0};
    const auto all_permutations2 = exponents2.allPermutations();
    const std::vector<GQCP::CartesianExponents> ref_all_permutations2 {GQCP::CartesianExponents(1, 0, 0), GQCP::CartesianExponents(0, 1, 0), GQCP::CartesianExponents(0, 0, 1)};
    BOOST_CHECK(all_permutations2 == ref_all_permutations2);


    // Example 3.
    const GQCP::CartesianExponents exponents3 {1, 1, 0};
    const auto all_permutations3 = exponents3.allPermutations();
    const std::vector<GQCP::CartesianExponents> ref_all_permutations3 {GQCP::CartesianExponents(1, 1, 0), GQCP::CartesianExponents(1, 0, 1), GQCP::CartesianExponents(0, 1, 1)};
    BOOST_CHECK(all_permutations3 == ref_all_permutations3);
}
