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
#define BOOST_TEST_MODULE "FrozenFockSpace"
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


#include "FockSpace/FrozenFockSpace.hpp"




BOOST_AUTO_TEST_CASE ( FrozenFockSpace_constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::FrozenFockSpace (10, 5, 1));
    BOOST_CHECK_NO_THROW(GQCP::FrozenFockSpace (10, 5, 2));
    BOOST_CHECK_NO_THROW(GQCP::FrozenFockSpace (10, 5, 3));
    BOOST_CHECK_NO_THROW(GQCP::FrozenFockSpace (10, 5, 4));
    BOOST_CHECK_NO_THROW(GQCP::FrozenFockSpace (10, 5, 5));
}

BOOST_AUTO_TEST_CASE ( FrozenFockSpace_dimension ) {

    BOOST_CHECK_EQUAL(GQCP::FrozenFockSpace::calculateDimension(10, 1), 10);
    BOOST_CHECK_EQUAL(GQCP::FrozenFockSpace::calculateDimension(6, 2), 15);
    BOOST_CHECK_EQUAL(GQCP::FrozenFockSpace::calculateDimension(8, 3), 56);
}


BOOST_AUTO_TEST_CASE ( vertex_weights_K5_N3 ) {

    // Let's test an addressing scheme for K=5 and N=3 (5 MOs and 3 alpha electrons)
    GQCP::FrozenFockSpace fock_space = GQCP::FrozenFockSpace(5, 3);

    GQCP::Matrixu ref_vertex_weights = {{1, 0, 0, 0},
                                        {1, 1, 0, 0},
                                        {1, 2, 1, 0},
                                        {0, 3, 3, 1},
                                        {0, 0, 6, 4},
                                        {0, 0, 0, 10}};
    BOOST_CHECK(ref_vertex_weights == fock_space.get_vertex_weights());
}


BOOST_AUTO_TEST_CASE ( iterateToNextUnoccupiedOrbital ) {

    GQCP::FrozenFockSpace fock_space (5, 3);
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

    GQCP::FrozenFockSpace fock_space (5, 3);
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


BOOST_AUTO_TEST_CASE ( shiftToPreviousOrbital_signed ) {

    GQCP::FrozenFockSpace fock_space (5, 3);
    GQCP::ONV onv = fock_space.makeONV(6);  // 10110

    size_t address_shift = 0;
    int sign = 1;

    // test shift if we plan on creating one electron and start from orbital index 2
    size_t e = 0;  // count starts at 0 (translates to orbital index 1)
    size_t q = 1;  // index starts at orbital index 1

    //  Index 1 is occupied an thus its weight shall shift if we create before an electron on a smaller index.
    //  In this instance the electron weight at index 1 should be shifted.
    //  Initial weight contribution is 1,
    //  This should be shifted to 0, the difference is 1.
    //  The sign changes (flips once)
    fock_space.shiftUntilPreviousUnoccupiedOrbital<1>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == -1);
    BOOST_CHECK(e == -1);
    BOOST_CHECK(q == 0);
    BOOST_CHECK(sign == -1);

    sign = 1;
    address_shift = 0;
    onv = fock_space.makeONV(9);  // 11100
    // test shift if we plan on creating two electrons and start from orbital index 2
    e = 0;  // count starts at 1 (translates to orbital index 2)
    q = 2;  // index starts at orbital index 2

    //  Index 2 is occupied an thus its weight shall shift if we create before an electron on a smaller index.
    //  In this instance the electron weight at index 2 should be shifted.
    //  Initial weight contribution is 2,
    //  This should be shifted to 0, the difference is 2.
    //  The sign changes (flips once)
    fock_space.shiftUntilPreviousUnoccupiedOrbital<2>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == -2);
    BOOST_CHECK(e == -1);
    BOOST_CHECK(q == 1);
    BOOST_CHECK(sign == -1);
}


BOOST_AUTO_TEST_CASE ( coupling_count ) {

    GQCP::FrozenFockSpace fock_space (5, 3);
    GQCP::ONV onv = fock_space.makeONV(3);  // 01110

    // We only count couplings with larger addresses

    BOOST_CHECK(fock_space.countOneElectronCouplings(onv) == 3); // 11100, 11010, 10110
    BOOST_CHECK(fock_space.countTwoElectronCouplings(onv) == 3+3);  // 11100, 11010, 10110, 11001, 10101, 10011

    onv = fock_space.makeONV(0);  // 00111

    BOOST_CHECK(fock_space.countOneElectronCouplings(onv) == 6);
    BOOST_CHECK(fock_space.countTwoElectronCouplings(onv) == 6+3); // all of them


    // test whether the total count matches that of individual counts of all ONVs in the Fock space.
    GQCP::FrozenFockSpace fock_space2 (16,8);

    size_t coupling_count1 = 0;
    size_t coupling_count2 = 0;
    onv = fock_space2.makeONV(0);  // spin string with address 0
    for (size_t I = 0; I < fock_space2.get_dimension(); I++) {  // I_alpha loops over all addresses of alpha spin strings
        if (I > 0) {
            fock_space2.setNextONV(onv);
        }
        coupling_count1 += fock_space2.countOneElectronCouplings(onv);
        coupling_count2 += fock_space2.countTwoElectronCouplings(onv);
    }

    BOOST_CHECK(2*coupling_count1 == fock_space2.countTotalOneElectronCouplings());
    BOOST_CHECK(2*coupling_count2 == fock_space2.countTotalTwoElectronCouplings());
}
