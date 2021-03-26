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

#define BOOST_TEST_MODULE "LeviCivitaTensor"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/LeviCivitaTensor.hpp"


/**
 *  Check if the Levi-Civita tensor's elements are correctly implemented.
 */
BOOST_AUTO_TEST_CASE(elements) {

    const GQCP::LeviCivitaTensor<double> epsilon {};

    BOOST_CHECK(epsilon(0, 0, 0) == 0);
    BOOST_CHECK(epsilon(0, 0, 1) == 0);
    BOOST_CHECK(epsilon(0, 0, 2) == 0);
    BOOST_CHECK(epsilon(0, 1, 0) == 0);
    BOOST_CHECK(epsilon(0, 1, 1) == 0);
    BOOST_CHECK(epsilon(0, 1, 2) == 1.0);
    BOOST_CHECK(epsilon(0, 2, 0) == 0);
    BOOST_CHECK(epsilon(0, 2, 1) == -1.0);
    BOOST_CHECK(epsilon(0, 2, 2) == 0);

    BOOST_CHECK(epsilon(1, 0, 0) == 0);
    BOOST_CHECK(epsilon(1, 0, 1) == 0);
    BOOST_CHECK(epsilon(1, 0, 2) == -1.0);
    BOOST_CHECK(epsilon(1, 1, 0) == 0);
    BOOST_CHECK(epsilon(1, 1, 1) == 0);
    BOOST_CHECK(epsilon(1, 1, 2) == 0);
    BOOST_CHECK(epsilon(1, 2, 0) == 1.0);
    BOOST_CHECK(epsilon(1, 2, 1) == 0);
    BOOST_CHECK(epsilon(1, 2, 2) == 0);

    BOOST_CHECK(epsilon(2, 0, 0) == 0);
    BOOST_CHECK(epsilon(2, 0, 1) == 1.0);
    BOOST_CHECK(epsilon(2, 0, 2) == 0);
    BOOST_CHECK(epsilon(2, 1, 0) == -1.0);
    BOOST_CHECK(epsilon(2, 1, 1) == 0);
    BOOST_CHECK(epsilon(2, 1, 2) == 0);
    BOOST_CHECK(epsilon(2, 2, 0) == 0);
    BOOST_CHECK(epsilon(2, 2, 1) == 0);
    BOOST_CHECK(epsilon(2, 2, 2) == 0);
}


/**
 *  Check if nonZeroIndex is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(nonZeroIndex) {

    const GQCP::LeviCivitaTensor<double> epsilon {};

    BOOST_CHECK(epsilon.nonZeroIndex(0, 1) == 2);
    BOOST_CHECK(epsilon.nonZeroIndex(1, 0) == 2);

    BOOST_CHECK(epsilon.nonZeroIndex(0, 2) == 1);
    BOOST_CHECK(epsilon.nonZeroIndex(2, 0) == 1);

    BOOST_CHECK(epsilon.nonZeroIndex(1, 2) == 0);
    BOOST_CHECK(epsilon.nonZeroIndex(1, 2) == 0);
}
