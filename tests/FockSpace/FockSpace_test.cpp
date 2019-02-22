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


BOOST_AUTO_TEST_CASE ( shiftToPreviousOrbital_signed ) {

    GQCP::FockSpace fock_space (5, 3);
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

    GQCP::FockSpace fock_space (5, 3);
    GQCP::ONV onv = fock_space.makeONV(3);  // 01110

    // We only count couplings with larger addresses

    BOOST_CHECK(fock_space.countOneElectronCouplings(onv) == 3); // 11100, 11010, 10110
    BOOST_CHECK(fock_space.countTwoElectronCouplings(onv) == 3+3);  // 11100, 11010, 10110, 11001, 10101, 10011

    onv = fock_space.makeONV(0);  // 00111

    BOOST_CHECK(fock_space.countOneElectronCouplings(onv) == 6);
    BOOST_CHECK(fock_space.countTwoElectronCouplings(onv) == 6+3); // all of them


    // test whether the total count matches that of individual counts of all ONVs in the Fock space.
    GQCP::FockSpace fock_space2 (16,8);

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

BOOST_AUTO_TEST_CASE ( ONV_address_setNext_fullspace ) {

    // Here we will test a full permutation through a Fock space of K = 15, N = 5
    GQCP::FockSpace fock_space (15, 5);

    // Retrieve the first ONV of the Fock space
    GQCP::ONV onv_test = fock_space.makeONV(0);

    const size_t dimension_fock_space = 3003;
    bool is_correct = true;  // variable that is updated to false if an unexpected result occurs

    // Iterate through the Fock space in reverse lexicographical order and test whether address matches
    for (size_t i = 0; i < dimension_fock_space; i++) {

        // Tests address
        if (i != fock_space.getAddress(onv_test)) {
            is_correct = false;
        }

        // transforms the given ONV to the next ONV in the Fock space
        if (i < dimension_fock_space - 1) {
            fock_space.setNextONV(onv_test);
        }
    }

    // Checks if no unexpected results occured in a full iteration
    BOOST_CHECK(is_correct);
}


BOOST_AUTO_TEST_CASE ( FockSpace_getAddress ) {

    GQCP::FockSpace fock_space (6, 3);

    // The address of the string "010011" (19) should be 4
    GQCP::ONV onv (6, 3, 19);

    BOOST_CHECK_EQUAL(fock_space.getAddress(onv), 4);
}


BOOST_AUTO_TEST_CASE ( FockSpace_setNext ) {

    GQCP::FockSpace fock_space (5, 3);
    // K = 5, N = 3 <-> "00111"
    GQCP::ONV onv = fock_space.makeONV(0);
    // The lexical permutations are: "00111" (7), "01011" (11), "01101" (13), "01110" (14), etc.

    // Check permutations one after the other

    fock_space.setNextONV(onv);  // "01011" (11)
    BOOST_CHECK_EQUAL(onv.get_unsigned_representation(), 11);
    GQCP::VectorXs x1 (3);
    x1 << 0, 1, 3;
    BOOST_CHECK(x1.isApprox(onv.get_occupation_indices()));

    fock_space.setNextONV(onv);  // "01101" (13)
    BOOST_CHECK_EQUAL(onv.get_unsigned_representation(), 13);
    GQCP::VectorXs x2 (3);
    x2 << 0, 2, 3;
    BOOST_CHECK(x2.isApprox(onv.get_occupation_indices()));

    fock_space.setNextONV(onv);  // "01110" (14)
    BOOST_CHECK_EQUAL(onv.get_unsigned_representation(), 14);
    GQCP::VectorXs x3 (3);
    x3 << 1, 2, 3;
    BOOST_CHECK(x3.isApprox(onv.get_occupation_indices()));
}
