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

#define BOOST_TEST_MODULE "JacobiRotation"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/JacobiRotation.hpp"


/**
 *  Check the constructor of JacobiRotation.
 */
BOOST_AUTO_TEST_CASE(JacobiRotation_constructor) {

    // Check if a correct constructor works.
    BOOST_CHECK_NO_THROW(GQCP::JacobiRotation(3, 1, 0.5));

    // Check if we can't construct when p < q.
    BOOST_CHECK_THROW(GQCP::JacobiRotation(1, 3, 0.5), std::invalid_argument);
}
