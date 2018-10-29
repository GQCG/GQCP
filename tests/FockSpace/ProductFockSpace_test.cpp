// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#define BOOST_TEST_MODULE "ProductFockSpace"


#include "FockSpace/ProductFockSpace.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( ProductFockSpace_constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::ProductFockSpace (10, 5, 5));
}


BOOST_AUTO_TEST_CASE ( ProductFockSpace_dimension) {

    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(10, 1, 1), 100);
    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(6, 2, 2), 225);
    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(8, 3, 3), 3136);

    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(10, 2, 0), 45);
    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(6, 3, 1), 120);
    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(8, 4, 2), 1960);
}
