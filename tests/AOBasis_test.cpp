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
#define BOOST_TEST_MODULE "AOBasis"


#include "AOBasis.hpp"

#include "Molecule.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( AOBasis_constructor ) {

    // Check if we can construct an AOBasis object
    auto water = GQCP::Molecule::Readxyz("../tests/data/h2o.xyz");
    GQCP::AOBasis basis (water, "STO-3G");
}


BOOST_AUTO_TEST_CASE ( number_of_basis_functions ) {

    // Check the number of basis functions in water
    auto water = GQCP::Molecule::Readxyz("../tests/data/h2o.xyz");
    GQCP::AOBasis basis (water, "STO-3G");

    BOOST_CHECK_EQUAL(basis.get_number_of_basis_functions(), 7);
}
