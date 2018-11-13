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
#define BOOST_TEST_MODULE "FockSpace"


#include "FockSpace/FockSpace.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( FockSpace_constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::FockSpace (10, 5));
}


BOOST_AUTO_TEST_CASE ( FockSpace_dimension) {

    BOOST_CHECK_EQUAL(GQCP::FockSpace::calculateDimension(10, 1), 10);
    BOOST_CHECK_EQUAL(GQCP::FockSpace::calculateDimension(6, 2), 15);
    BOOST_CHECK_EQUAL(GQCP::FockSpace::calculateDimension(8, 3), 56);
}


BOOST_AUTO_TEST_CASE ( vertex_weights_K5_N3 ) {

    // Let's test an addressing scheme for K=5 and N=3 (5 MOs and 3 alpha electrons)
    GQCP::FockSpace fock_space = GQCP::FockSpace(5, 3);

    GQCP::Matrixu ref_vertex_weights = {{1, 0, 0, 0},
                                        {1, 1, 0, 0},
                                        {1, 2, 1, 0},
                                        {0, 3, 3, 1},
                                        {0, 0, 6, 4},
                                        {0, 0, 0, 10}};
    BOOST_CHECK(ref_vertex_weights == fock_space.get_vertex_weights());
}


BOOST_AUTO_TEST_CASE ( iterateToNextUnoccupiedOrbital ) {

    GQCP::FockSpace fock_space (5, 3);
    GQCP::ONV onv = fock_space.get_ONV(3);  // 01110

    size_t e = 0;  // count starts at 0 (one electron = 0, two electrons = 1, etc)
    size_t q = 2;  // index starts at the second occupied orbital (as if the first was annihilated)

    //  In this instance electron weights at index 2 and 3 should be shifted
    //  Initial weight contributions were 1 and 1 respectively.
    //  These should be shifted to 2 and 3 respectively, the difference is 1 and 2 respectively
    //  The total shift is thus 3
    size_t address_shift = fock_space.shiftAddressTillNextUnoccupiedOrbital(onv, q, e);

    BOOST_CHECK(address_shift == 3);
}
