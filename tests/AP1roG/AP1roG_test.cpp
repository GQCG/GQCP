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
#define BOOST_TEST_MODULE "AP1roG"


#include "AP1roG/AP1roG.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( default_constructor ) {
    GQCP::AP1roG ap1rog;
}


BOOST_AUTO_TEST_CASE ( constructor ) {

    GQCP::AP1roGGeminalCoefficients g (4, 6);
    GQCP::AP1roG ap1rog (g, 0.0);
}


BOOST_AUTO_TEST_CASE ( get_geminal_coefficients ) {

    GQCP::AP1roGGeminalCoefficients g (4, 6);
    GQCP::AP1roG ap1rog (g, 0.0);
    ap1rog.get_geminal_coefficients();
}
