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
#define BOOST_TEST_MODULE "JacobiRotationParameters"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "JacobiRotationParameters.hpp"


BOOST_AUTO_TEST_CASE ( JacobiRotationParameters_constructor ) {

    // Check if a correct constructor works
    GQCP::JacobiRotationParameters jacobi_parameters (3, 1, 0.5);

    // Check if we can't construct when p < q
    BOOST_CHECK_THROW(GQCP::JacobiRotationParameters(1, 3, 0.5), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( operator_stream ) {

    GQCP::JacobiRotationParameters jacobi_parameters (3, 1, 0.5);
    std::cout << jacobi_parameters << std::endl;
}
