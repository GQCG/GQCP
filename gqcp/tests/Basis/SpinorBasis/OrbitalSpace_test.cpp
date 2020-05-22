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


// Create some shortcuts to be used in the following tests.
const auto occ = GQCP::OccupationType::k_occupied;
const auto act = GQCP::OccupationType::k_active;
const auto virt = GQCP::OccupationType::k_virtual;


/**
 *  Check if the constructor and basic public methods work as expected.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    const std::vector<size_t> occupied_indices {0, 1, 2};
    const std::vector<size_t> active_indices {3};
    const std::vector<size_t> virtual_indices {4, 5, 6, 7};
    const std::vector<size_t> all_indices {0, 1, 2, 3, 4, 5, 6, 7};

    const GQCP::OrbitalSpace orbital_space {occupied_indices, active_indices, virtual_indices};

    BOOST_CHECK(orbital_space.indices() == all_indices);
    BOOST_CHECK(orbital_space.indices(occ) == occupied_indices);
    BOOST_CHECK(orbital_space.indices(act) == active_indices);
    BOOST_CHECK(orbital_space.indices(virt) == virtual_indices);

    BOOST_CHECK(orbital_space.numberOfOrbitals() == 8);
    BOOST_CHECK(orbital_space.numberOfOrbitals(occ) == 3);
    BOOST_CHECK(orbital_space.numberOfOrbitals(act) == 1);
    BOOST_CHECK(orbital_space.numberOfOrbitals(virt) == 4);
}


/**
 *  Check if the named constructor ::Implicit works as expected.
 */
BOOST_AUTO_TEST_CASE(Implicit) {


    // Check an occupied-only constructor.
    const std::vector<size_t> ref_occupied_indices1 {0, 1, 2, 3};

    const auto orbital_space1 = GQCP::OrbitalSpace::Implicit({{occ, 4}});
    BOOST_CHECK(orbital_space1.indices(occ) == ref_occupied_indices1);


    // Check an occupied-virtual constructor.
    const std::vector<size_t> ref_occupied_indices2 {0, 1, 2, 3};
    const std::vector<size_t> ref_virtual_indices2 {4, 5, 6};

    const auto orbital_space2 = GQCP::OrbitalSpace::Implicit({{occ, 4}, {virt, 3}});
    BOOST_CHECK(orbital_space2.indices(occ) == ref_occupied_indices2);
    BOOST_CHECK(orbital_space2.indices(virt) == ref_virtual_indices2);


    // Check an occupied-active-virtual constructor.
    const std::vector<size_t> ref_occupied_indices3 {0, 1, 2, 3};
    const std::vector<size_t> ref_active_indices3 {4, 5};
    const std::vector<size_t> ref_virtual_indices3 {6, 7, 8, 9, 10};


    const auto orbital_space3 = GQCP::OrbitalSpace::Implicit({{occ, 4}, {virt, 5}, {act, 2}});
    BOOST_CHECK(orbital_space3.indices(occ) == ref_occupied_indices3);
    BOOST_CHECK(orbital_space3.indices(act) == ref_active_indices3);
    BOOST_CHECK(orbital_space3.indices(virt) == ref_virtual_indices3);
}


/**
 *  Check if the isIndex() method works as expected.
 */
BOOST_AUTO_TEST_CASE(isIndex) {

    const std::vector<size_t> occupied_indices {0, 1, 2};
    const std::vector<size_t> active_indices {3};
    const std::vector<size_t> virtual_indices {4, 5, 6, 7};

    const GQCP::OrbitalSpace orbital_space {occupied_indices, active_indices, virtual_indices};


    BOOST_CHECK(orbital_space.isIndex(act, 3));

    BOOST_CHECK(orbital_space.isIndex(occ, 0));
    BOOST_CHECK(orbital_space.isIndex(occ, 1));
    BOOST_CHECK(orbital_space.isIndex(occ, 2));

    BOOST_CHECK(orbital_space.isIndex(virt, 4));
    BOOST_CHECK(orbital_space.isIndex(virt, 5));
    BOOST_CHECK(orbital_space.isIndex(virt, 6));
    BOOST_CHECK(orbital_space.isIndex(virt, 7));
}
