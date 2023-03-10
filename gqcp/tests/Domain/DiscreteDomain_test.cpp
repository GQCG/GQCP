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

#define BOOST_TEST_MODULE "DiscreteDomain"

#include <boost/test/unit_test.hpp>

#include "Domain/DiscreteDomain.hpp"
#include "Domain/SimpleDomain.hpp"


/**
 *  Test if the bitstring representation of the domain is returned correctly.
 */
BOOST_AUTO_TEST_CASE(bitstring_representation) {
    // Construct a discrete domain "" of size 8 where the indices 0, 2, 4 and 6 are occupied.
    // Note that the bitstring representation of the domain should be read from right to left.
    GQCP::VectorX<size_t> element_indices {0};
    element_indices << 0, 2, 4, 6;
    const GQCP::DiscreteDomain discrete_domain {element_indices, 8};
    BOOST_CHECK_EQUAL(discrete_domain.asString(), std::string("01010101"));
}


/**
 *  Test if the `addElement' functionality correctly adds an element to the discrete domain.
 */
BOOST_AUTO_TEST_CASE(add_element) {
    GQCP::VectorX<size_t> element_indices {0};
    GQCP::DiscreteDomain discrete_domain {element_indices, 8};

    // Add domain element at index 0.
    discrete_domain.addElement(0);
    element_indices << 0;
    BOOST_CHECK_EQUAL(discrete_domain(0), 1);
    BOOST_CHECK_EQUAL(discrete_domain.asString(), std::string("00000001"));
    BOOST_CHECK(discrete_domain.domainIndices().isApprox(element_indices, 1e-06));

    // Add domain element at index 3.
    discrete_domain.addElement(3);
    element_indices << 3;
    BOOST_CHECK_EQUAL(discrete_domain(3), 1);
    BOOST_CHECK_EQUAL(discrete_domain.asString(), std::string("00001001"));
    BOOST_CHECK(discrete_domain.domainIndices().isApprox(element_indices, 1e-06));

    // Check whether the number of domain elements is now 2.
    BOOST_CHECK_EQUAL(discrete_domain.numberOfElements(), 2);
}


/**
 *  Test if the `removeElement' functionality correctly removes an element from the discrete domain.
 */
BOOST_AUTO_TEST_CASE(remove_element) {
    GQCP::VectorX<size_t> element_indices {0};
    element_indices << 0, 5;
    GQCP::DiscreteDomain discrete_domain {element_indices, 8};

    // Add domain element at index 0.
    discrete_domain.removeElement(0);
    element_indices << 0;
    BOOST_CHECK_EQUAL(discrete_domain(0), 0);
    BOOST_CHECK_EQUAL(discrete_domain.asString(), std::string("00100000"));
    BOOST_CHECK(discrete_domain.domainIndices().isApprox(element_indices, 1e-06));

    // Add domain element at index 3.
    discrete_domain.addElement(3);
    element_indices << 3;
    BOOST_CHECK_EQUAL(discrete_domain(3), 1);
    BOOST_CHECK_EQUAL(discrete_domain.asString(), std::string("00000000"));
    BOOST_CHECK(discrete_domain.domainIndices().isApprox(element_indices, 1e-06));

    // Check whether the number of domain elements is now 2.
    BOOST_CHECK_EQUAL(discrete_domain.numberOfElements(), 0);
}
