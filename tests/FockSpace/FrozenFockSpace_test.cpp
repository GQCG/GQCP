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


BOOST_AUTO_TEST_CASE ( coupling_count ) {

    GQCP::FrozenFockSpace fock_space (7, 5, 2);
    GQCP::ONV onv = fock_space.makeONV(3);

    // We only count couplings with larger addresses

    BOOST_CHECK(fock_space.countOneElectronCouplings(onv) == 3);
    BOOST_CHECK(fock_space.countTwoElectronCouplings(onv) == 3+3);

    onv = fock_space.makeONV(0);

    BOOST_CHECK(fock_space.countOneElectronCouplings(onv) == 6);
    BOOST_CHECK(fock_space.countTwoElectronCouplings(onv) == 6+3);


    // test whether the total count matches that of individual counts of all ONVs in the Fock space.
    GQCP::FrozenFockSpace fock_space2 (18, 10, 2);

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


BOOST_AUTO_TEST_CASE ( ulongNextPermutation_getAddress_calculateRepresentation ) {

    size_t K = 5;
    size_t N = 3;
    size_t X = 1;  // number of frozen orbitals/electrons

    GQCP::FrozenFockSpace frozen_fock_space(K, N, X);

    // "00111", "01011", "01101", "10011", "10101", "11001"
    size_t bit_frozen_fock_space[6] = {7, 11, 13, 19, 21, 25};

    size_t bit_onv = bit_frozen_fock_space[0];
    BOOST_CHECK(frozen_fock_space.getAddress(bit_onv) == 0);
    for (size_t i = 1; i < frozen_fock_space.get_dimension(); i++) {
        bit_onv = frozen_fock_space.ulongNextPermutation(bit_onv);
        BOOST_CHECK(bit_onv == bit_frozen_fock_space[i]);
        BOOST_CHECK(frozen_fock_space.getAddress(bit_onv) == i);
        BOOST_CHECK(frozen_fock_space.calculateRepresentation(i) == bit_onv);
    }
}


BOOST_AUTO_TEST_CASE ( address_setNext_frozen_space ) {

    // Here we will test a set of permutations through a frozen Fock space of K = 15, N = 5, X = 2
    GQCP::FrozenFockSpace fock_space (15, 5, 2);

    // Retrieve the first ONV of the Fock space
    GQCP::ONV onv_test = fock_space.makeONV(0);

    const size_t permutations = 285;
    bool is_correct = true;  // variable that is updated to false if an unexpected result occurs

    // Iterate through the Fock space in reverse lexicographical order and test whether address matches
    for (size_t i = 0; i < permutations; i++) {

        // Tests address
        if (i != fock_space.getAddress(onv_test)) {
            is_correct = false;
        }

        // transforms the given ONV to the next ONV in the Fock space
        if (i < permutations - 1) {
            fock_space.setNextONV(onv_test);
        }
    }

    // Checks if no unexpected results occurred in a full iteration
    BOOST_CHECK(is_correct);
}
