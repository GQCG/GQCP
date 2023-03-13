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
    const std::vector<size_t> element_indices = {0, 2, 4, 6};

    GQCP::DiscreteDomain discrete_domain_from_indices(element_indices, 8);
    BOOST_CHECK_EQUAL(discrete_domain_from_indices.asString(), std::string("10101010"));

    // The same as above, but with passing an unsigned integer to the constructor.
    // GQCP::DiscreteDomain discrete_domain_from_unsigned_int(85, 8);
    // BOOST_CHECK_EQUAL(discrete_domain_from_unsigned_int.asString(), std::string("10101010"));
}


/**
 *  Test if the `addElement' functionality correctly adds an element to the discrete domain.
 */
BOOST_AUTO_TEST_CASE(add_element) {
    // Start from an empty domain.
    GQCP::DiscreteDomain discrete_domain {0, 8};

    // Add domain element at index 0.
    discrete_domain.addElement(0);
    std::vector<size_t> element_indices {0};
    BOOST_CHECK_EQUAL(discrete_domain(0), 1);
    BOOST_CHECK_EQUAL(discrete_domain.asString(), std::string("10000000"));
    BOOST_CHECK_EQUAL_COLLECTIONS(discrete_domain.domainIndices().begin(), discrete_domain.domainIndices().end(), element_indices.begin(), element_indices.end());

    // Add domain element at index 3.
    discrete_domain.addElement(3);
    element_indices.push_back(3);
    BOOST_CHECK_EQUAL(discrete_domain(3), 1);
    BOOST_CHECK_EQUAL(discrete_domain.asString(), std::string("10010000"));
    BOOST_CHECK_EQUAL_COLLECTIONS(discrete_domain.domainIndices().begin(), discrete_domain.domainIndices().end(), element_indices.begin(), element_indices.end());

    // Check whether the number of domain elements is now 2.
    BOOST_CHECK_EQUAL(discrete_domain.numberOfElements(), 2);
}


/**
 *  Test if the `removeElement' functionality correctly removes an element from the discrete domain.
 */
BOOST_AUTO_TEST_CASE(remove_element) {
    std::vector<size_t> element_indices {0, 5};
    GQCP::DiscreteDomain discrete_domain {element_indices, 8};

    // Remove domain element at index 5.
    discrete_domain.removeElement(5);
    element_indices.pop_back();
    BOOST_CHECK_EQUAL(discrete_domain(5), 0);
    BOOST_CHECK_EQUAL(discrete_domain.asString(), std::string("10000000"));
    BOOST_CHECK_EQUAL_COLLECTIONS(discrete_domain.domainIndices().begin(), discrete_domain.domainIndices().end(), element_indices.begin(), element_indices.end());

    // Remove domain element at index 0.
    discrete_domain.removeElement(0);
    element_indices.pop_back();
    BOOST_CHECK_EQUAL(discrete_domain(0), 0);
    BOOST_CHECK_EQUAL(discrete_domain.asString(), std::string("00000000"));
    BOOST_CHECK_EQUAL_COLLECTIONS(discrete_domain.domainIndices().begin(), discrete_domain.domainIndices().end(), element_indices.begin(), element_indices.end());

    // Check whether the number of domain elements is now 2.
    BOOST_CHECK_EQUAL(discrete_domain.numberOfElements(), 0);
}


/**
 *  Test whether the discrete domain assigns the correct unsigned representation associated with the given domain indices and bitstring.
 */
BOOST_AUTO_TEST_CASE(unsigned_representation) {

    // Discrete domain with bitstring representation "01001010".
    const std::vector<size_t> element_indices1 {1, 3, 6};
    GQCP::DiscreteDomain discrete_domain {element_indices1, 8};
    BOOST_CHECK_EQUAL(discrete_domain.unsignedRepresentation(), 74);

    const std::vector<size_t> element_indices2 {0, 1, 2, 3, 4, 5, 6, 7};
    discrete_domain = GQCP::DiscreteDomain {element_indices2, 8};
    BOOST_CHECK_EQUAL(discrete_domain.unsignedRepresentation(), 255);
}
