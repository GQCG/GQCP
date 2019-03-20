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

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "Basis/CartesianExponents.hpp"


BOOST_AUTO_TEST_CASE ( angular_momentum ) {

    auto exps = GQCP::CartesianExponents({2, 3, 4});  // exponents

    BOOST_CHECK_EQUAL(exps.angularMomentum(), 9);
}


BOOST_AUTO_TEST_CASE ( operator_lt ) {

    auto exps1 = GQCP::CartesianExponents({1, 0, 0});  // exponents
    auto exps2 = GQCP::CartesianExponents({2, 0, 0});  // exponents
    auto exps3 = GQCP::CartesianExponents({1, 1, 0});  // exponents
    auto exps4 = GQCP::CartesianExponents({0, 1, 1});  // exponents

    BOOST_CHECK(exps1 < exps2);
    BOOST_CHECK(exps1 < exps3);
    BOOST_CHECK(exps1 < exps4);

    BOOST_CHECK(exps2 < exps3);
    BOOST_CHECK(exps2 < exps4);
    BOOST_CHECK(exps3 < exps4);
}
