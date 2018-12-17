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
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


#include "FockSpace/FockSpace.hpp"




BOOST_AUTO_TEST_CASE ( FockSpace_constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::FockSpace (10, 5));
}


BOOST_AUTO_TEST_CASE ( expansions ) {

    // Check the BaseFockSpace expansion functions
    GQCP::FockSpace fock_space (8, 3);

    Eigen::VectorXd hartree_fock_expansion = fock_space.HartreeFockExpansion();
    BOOST_CHECK(std::abs(hartree_fock_expansion.norm() - 1.0) < 1.0e-12);  // check if normalized
    BOOST_CHECK(std::abs(hartree_fock_expansion(0) - 1.0) < 1.0e-12);  // the Hartree-Fock determinant should be the first one

    Eigen::VectorXd random_expansion = fock_space.randomExpansion();
    BOOST_CHECK(std::abs(random_expansion.norm() - 1.0) < 1.0e-12);  // check if normalized

    Eigen::VectorXd constant_expansion = fock_space.constantExpansion();
    BOOST_CHECK(std::abs(constant_expansion.norm() - 1.0) < 1.0e-12);  // check if normalized
}


BOOST_AUTO_TEST_CASE ( FockSpace_dimension ) {

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
    GQCP::ONV onv = fock_space.makeONV(3);  // 01110

    size_t address_shift = 0;
    // test shift if we annihilate one electron and start from orbital index 2
    size_t e = 1;  // count starts at 1 (translates to orbital index 2)
    size_t q = 2;  // index starts at orbital index 2

    //  In this instance electron weights at index 2 and 3 should be shifted.
    //  Initial weight contributions were 1 and 1 respectively,
    //  these should be shifted to 2 and 3 respectively, the difference is 1 and 2 respectively.
    //  The total shift is thus 3
    fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address_shift, q, e);

    BOOST_CHECK(address_shift == 3);
    BOOST_CHECK(e == 3);
    BOOST_CHECK(q == 4);

    // test shift if we annihilate two electrons and start from orbital index 3
    e = 2;  // count starts at 2 (translates to orbital index 3)
    q = 3;  // index starts at orbital index 3

    //  In this instance electron weights at index 3 should be shifted.
    //  The initial weight contribution was 1,
    //  this should be shifted to 3, the difference is 2
    //  The total shift is thus 2
    address_shift = 0;
    fock_space.shiftUntilNextUnoccupiedOrbital<2>(onv, address_shift, q, e);

    BOOST_CHECK(address_shift == 2);
    BOOST_CHECK(e == 3);
    BOOST_CHECK(q == 4);
}


BOOST_AUTO_TEST_CASE ( iterateToNextUnoccupiedOrbital_signed ) {

    GQCP::FockSpace fock_space (5, 3);
    GQCP::ONV onv = fock_space.makeONV(3);  // 01110

    size_t address_shift = 0;
    int sign = 1;
    // test shift if we annihilate one electron and start from orbital index 2
    size_t e = 1;  // count starts at 1 (translates to orbital index 2)
    size_t q = 2;  // index starts at orbital index 2

    //  In this instance electron weights at index 2 and 3 should be shifted.
    //  Initial weight contributions were 1 and 1 respectively,
    //  these should be shifted to 2 and 3 respectively, the difference is 1 and 2 respectively.
    //  The total shift is thus 3, and the sign should remain the same (flips twice)
    fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == 3);
    BOOST_CHECK(e == 3);
    BOOST_CHECK(q == 4);
    BOOST_CHECK(sign == 1);

    // test shift if we annihilate two electrons and start from orbital index 3
    e = 2;  // count starts at 2 (translates to orbital index 3)
    q = 3;  // index starts at orbital index 3

    //  In this instance electron weights at index 3 should be shifted.
    //  The initial weight contribution was 1,
    //  this should be shifted to 3, the difference is 2
    //  The total shift is thus 2, sign should flip once
    address_shift = 0;
    fock_space.shiftUntilNextUnoccupiedOrbital<2>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == 2);
    BOOST_CHECK(e == 3);
    BOOST_CHECK(q == 4);
    BOOST_CHECK(sign == -1);
}
