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
#define BOOST_TEST_MODULE "CartesianExponents"

#include <boost/test/included/unit_test.hpp>

#include "Basis/CartesianExponents.hpp"


BOOST_AUTO_TEST_CASE ( angular_momentum ) {

    GQCP::CartesianExponents exps (2, 3, 4);  // exps: exponents

    BOOST_CHECK_EQUAL(exps.angularMomentum(), 9);
}


BOOST_AUTO_TEST_CASE ( operator_smaller ) {

    GQCP::CartesianExponents exps1 (1, 0, 0);  // exps: exponents
    GQCP::CartesianExponents exps2 (2, 0, 0);  // exps: exponents
    GQCP::CartesianExponents exps3 (1, 1, 0);  // exps: exponents
    GQCP::CartesianExponents exps4 (0, 1, 1);  // exps: exponents

    BOOST_CHECK(exps1 < exps2);
    BOOST_CHECK(exps1 < exps3);
    BOOST_CHECK(exps1 < exps4);

    BOOST_CHECK(exps2 < exps3);
    BOOST_CHECK(exps2 < exps4);
    BOOST_CHECK(exps3 < exps4);
}


BOOST_AUTO_TEST_CASE ( operator_equals ) {

    GQCP::CartesianExponents exps1 (1, 0, 0);  // exps: exponents
    GQCP::CartesianExponents exps2 (1, 0, 0);  // exps: exponents
    GQCP::CartesianExponents exps3 (2, 0, 0);  // exps: exponents

    BOOST_CHECK(exps1 == exps2);
    BOOST_CHECK(!(exps1 == exps3));
}


BOOST_AUTO_TEST_CASE ( allPermutations ) {

    GQCP::CartesianExponents exps1 (2, 3, 4);  // exps: exponents
    auto all_permutations1 = exps1.allPermutations();
    std::vector<GQCP::CartesianExponents> ref_all_permutations1 {
        GQCP::CartesianExponents(4, 3, 2), GQCP::CartesianExponents(4, 2, 3),
        GQCP::CartesianExponents(3, 4, 2), GQCP::CartesianExponents(3, 2, 4),
        GQCP::CartesianExponents(2, 4, 3), GQCP::CartesianExponents(2, 3, 4)
    };
    BOOST_CHECK(all_permutations1 == ref_all_permutations1);


    GQCP::CartesianExponents exps2 (1, 0, 0);  // exps: exponents
    auto all_permutations2 = exps2.allPermutations();
    std::vector<GQCP::CartesianExponents> ref_all_permutations2 {GQCP::CartesianExponents(1, 0, 0), GQCP::CartesianExponents(0, 1, 0), GQCP::CartesianExponents(0, 0, 1)};
    BOOST_CHECK(all_permutations2 == ref_all_permutations2);


    GQCP::CartesianExponents exps3 (1, 1, 0);  // exps: exponents
    auto all_permutations3 = exps3.allPermutations();
    std::vector<GQCP::CartesianExponents> ref_all_permutations3 {GQCP::CartesianExponents(1, 1, 0), GQCP::CartesianExponents(1, 0, 1), GQCP::CartesianExponents(0, 1, 1)};
    BOOST_CHECK(all_permutations3 == ref_all_permutations3);
}
