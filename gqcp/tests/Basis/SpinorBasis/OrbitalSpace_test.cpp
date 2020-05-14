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

#define BOOST_TEST_MODULE "OrbitalSpace"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/OrbitalSpace.hpp"


/**
 *  Check if the constructor and basic public methods work as expected.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    const std::vector<size_t> occupied_indices {0, 1, 2};
    const std::vector<size_t> active_indices {3};
    const std::vector<size_t> virtual_indices {4, 5, 6, 7};
    const std::vector<size_t> all_indices {0, 1, 2, 3, 4, 5, 6, 7};

    const GQCP::OrbitalSpace orbital_space {occupied_indices, active_indices, virtual_indices};

    BOOST_CHECK(orbital_space.activeIndices() == active_indices);
    BOOST_CHECK(orbital_space.allIndices() == all_indices);
    BOOST_CHECK(orbital_space.occupiedIndices() == occupied_indices);
    BOOST_CHECK(orbital_space.virtualIndices() == virtual_indices);

    BOOST_CHECK(orbital_space.numberOfActiveOrbitals() == 1);
    BOOST_CHECK(orbital_space.numberOfOccupiedOrbitals() == 3);
    BOOST_CHECK(orbital_space.numberOfOrbitals() == 8);
    BOOST_CHECK(orbital_space.numberOfVirtualOrbitals() == 4);
}


/**
 *  Check if the named constructor ::Occupied works as expected
 */
BOOST_AUTO_TEST_CASE(Occupied) {

    const std::vector<size_t> ref_occupied_indices {0, 1, 2, 3};

    BOOST_CHECK(GQCP::OrbitalSpace::Occupied(4).occupiedIndices() == ref_occupied_indices);
}


/**
 *  Check if the named constructor ::OccupiedVirtual works as expected
 */
BOOST_AUTO_TEST_CASE(OccupiedVirtual) {

    const std::vector<size_t> ref_occupied_indices {0, 1, 2, 3};
    const std::vector<size_t> ref_virtual_indices {4, 5, 6};

    const auto orbital_space = GQCP::OrbitalSpace::OccupiedVirtual(4, 7);
    BOOST_CHECK(orbital_space.occupiedIndices() == ref_occupied_indices);
    BOOST_CHECK(orbital_space.virtualIndices() == ref_virtual_indices);
}


/**
 *  Check if the isIndex<...> methods work as expected.
 */
BOOST_AUTO_TEST_CASE(isIndex) {

    const std::vector<size_t> occupied_indices {0, 1, 2};
    const std::vector<size_t> active_indices {3};
    const std::vector<size_t> virtual_indices {4, 5, 6, 7};

    const GQCP::OrbitalSpace orbital_space {occupied_indices, active_indices, virtual_indices};


    BOOST_CHECK(orbital_space.isIndexActive(3));

    BOOST_CHECK(orbital_space.isIndexOccupied(0));
    BOOST_CHECK(orbital_space.isIndexOccupied(1));
    BOOST_CHECK(orbital_space.isIndexOccupied(2));

    BOOST_CHECK(orbital_space.isIndexVirtual(4));
    BOOST_CHECK(orbital_space.isIndexVirtual(5));
    BOOST_CHECK(orbital_space.isIndexVirtual(6));
    BOOST_CHECK(orbital_space.isIndexVirtual(7));
}
