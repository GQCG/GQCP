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
 *  Test if the expectation value of a one-electron operator in different orbital bases is the same.
 */
BOOST_AUTO_TEST_CASE(bitstring_representation) {

    // Construct a discrete domain "" of size 8 where the indices 0, 2, 4 and 6 are occupied.
    // Note that the bitstring representation of the domain should be read from right to left.
    GQCP::VectorX<size_t> element_indices {0};
    element_indices << 0, 2, 4, 6;
    const size_t M = 8;
    const GQCP::DiscreteDomain discrete_domain {element_indices, M};
    std::cout << discrete_domain.asBitstring();
}
