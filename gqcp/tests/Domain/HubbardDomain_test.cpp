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

#define BOOST_TEST_MODULE "HubbardDomain"

#include <boost/test/unit_test.hpp>

#include "Domain/HubbardDomain.hpp"


/**
 *  Test whether the overlap between the Hubbard domain and a spin-unresolved ONV is calculated correctly.
 */
BOOST_AUTO_TEST_CASE(overlap_with_spinunresolved_onv) {
    // Create a Hubbard domain with bitstring representation 00001011 (read right to left).
    const GQCP::HubbardDomain domain {11, 8};
    // Create a spin-unresolved ONV with bitstring representation 00001011 (read right to left).
    GQCP::SpinUnresolvedONV onv {8, 3, 11};
    // The number of overlapping bits should be 3.
    BOOST_CHECK_EQUAL(domain.overlapWith(onv), 3);

    // Remove an electron from the ONV at position 1, the bitstring representation of the ONV is now 00001001 (read right to left).
    const auto annihilated = onv.annihilate(1);
    // The number of overlapping bits should now be 2.
    BOOST_CHECK_EQUAL(domain.overlapWith(onv), 2);

    // Create an electron at position 7, the bitstring representation of the ONV is now 10001001 (read right to left).
    const auto created = onv.create(7);
    // The number of overlapping bits should still be 2 since the domain does not contain the element at position 7.
    BOOST_CHECK_EQUAL(domain.overlapWith(onv), 2);

    // Create an empty Hubbard domain with bitstring representation 00000000 (read right to left).
    const GQCP::HubbardDomain empty_domain {0, 8};
    // Create a spin-unresolved ONV with bitstring representation 11111111 (read right to left).
    GQCP::SpinUnresolvedONV full_onv {8, 8, 255};
    // There cannot be an overlap when the domain is empty.
    BOOST_CHECK_EQUAL(empty_domain.overlapWith(full_onv), 0);
}


/**
 *  Test whether the overlap between the Hubbard domain and a spin-resolved ONV is calculated correctly.
 */
BOOST_AUTO_TEST_CASE(overlap_with_spinresolved_onv) {
    // Create a Hubbard domain with bitstring representation 00001011 (read right to left).
    const GQCP::HubbardDomain domain {11, 8};
    // Create a spin-resolved ONV where the alpha and beta part have bitstring representation 00101101 and 00010010 respectively (read right to left).
    const GQCP::SpinResolvedONV onv {GQCP::SpinUnresolvedONV {8, 4, 45}, GQCP::SpinUnresolvedONV {8, 2, 18}};
    // The number of overlapping bits should be 2 and 1 for both the alpha and beta part of the ONV respectively.
    BOOST_CHECK_EQUAL(domain.overlapWith(onv).alpha(), 2);
    BOOST_CHECK_EQUAL(domain.overlapWith(onv).beta(), 1);
}
