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

#define BOOST_TEST_MODULE "SpinResolvedOrbitalSpace"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Basis/SpinorBasis/SpinResolvedOrbitalSpace.hpp"


// Create some shortcuts to be used in the following tests.
const auto occ = GQCP::OccupationType::k_occupied;
const auto act = GQCP::OccupationType::k_active;
const auto virt = GQCP::OccupationType::k_virtual;

const auto alpha = GQCP::Spin::alpha;
const auto beta = GQCP::Spin::beta;


/**
 *  Check whether the constructor works as expected.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    // Create potential sets of indices for all occupation types.
    const std::vector<size_t> occupied_indices {0, 1, 2};
    const std::vector<size_t> active_indices {3};
    const std::vector<size_t> virtual_indices {4, 5, 6, 7};
    const std::vector<size_t> all_indices {0, 1, 2, 3, 4, 5, 6, 7};

    // Construct an alpha and beta orbital space.
    const GQCP::OrbitalSpace orbital_space_a {occupied_indices, active_indices, virtual_indices};
    const GQCP::OrbitalSpace orbital_space_b {occupied_indices, active_indices, virtual_indices};

    //Create the spin resolved orbital space from the alpha and beta component.
    const GQCP::SpinResolvedOrbitalSpace orbital_space {orbital_space_a, orbital_space_b};

    // Test whether or not the indices match. Some are tested for alpha, some for beta.
    BOOST_CHECK(orbital_space.alpha().indices() == all_indices);
    BOOST_CHECK(orbital_space.beta().indices(occ) == occupied_indices);
    BOOST_CHECK(orbital_space.alpha().indices(act) == active_indices);
    BOOST_CHECK(orbital_space.beta().indices(virt) == virtual_indices);

    BOOST_CHECK(orbital_space.alpha().numberOfOrbitals() == 8);
    BOOST_CHECK(orbital_space.beta().numberOfOrbitals(occ) == 3);
    BOOST_CHECK(orbital_space.alpha().numberOfOrbitals(act) == 1);
    BOOST_CHECK(orbital_space.beta().numberOfOrbitals(virt) == 4);
}


/**
 *  Check whether the initializing a mixed representable matrix object works as expected.
 */
BOOST_AUTO_TEST_CASE(mixed_matrix) {

    // Create occupied and virtual test indices for the alpha component.
    const std::vector<size_t> occupied_indices_a {0, 1, 2, 3};
    const std::vector<size_t> virtual_indices_a {4, 5, 6, 7, 8};

    // Do the same for the beta component.
    const std::vector<size_t> occupied_indices_b {0, 1, 2};
    const std::vector<size_t> virtual_indices_b {3, 4, 5};

    // Use the indices to create the two components and the total spin resolved orbital space.
    const GQCP::OrbitalSpace orbital_space_a {occupied_indices_a, virtual_indices_a};
    const GQCP::OrbitalSpace orbital_space_b {occupied_indices_b, virtual_indices_b};
    const GQCP::SpinResolvedOrbitalSpace orbital_space {orbital_space_a, orbital_space_b};

    // Zero-initialize a mixed spin matrix object.
    const auto M = orbital_space.initializeMixedRepresentableObjectFor<double>(occ, alpha, virt, beta);

    // Test whether the dimensions of the mixed object match the correct number of indices of the chosen occupation type and spin.
    BOOST_CHECK(M.asMatrix().rows() == occupied_indices_a.size());
    BOOST_CHECK(M.asMatrix().cols() == virtual_indices_b.size());
}


/**
 *  Check whether the initializing a mixed representable tensor object works as expected.
 */
BOOST_AUTO_TEST_CASE(mixed_tensor) {

    // Create occupied and virtual test indices for the alpha component.
    const std::vector<size_t> occupied_indices_a {0, 1, 2, 3};
    const std::vector<size_t> virtual_indices_a {4, 5, 6, 7, 8};

    // Do the same for the beta component.
    const std::vector<size_t> occupied_indices_b {0, 1, 2};
    const std::vector<size_t> virtual_indices_b {3, 4, 5};

    // Use the indices to create the two components and the total spin resolved orbital space.
    const GQCP::OrbitalSpace orbital_space_a {occupied_indices_a, virtual_indices_a};
    const GQCP::OrbitalSpace orbital_space_b {occupied_indices_b, virtual_indices_b};
    const GQCP::SpinResolvedOrbitalSpace orbital_space {orbital_space_a, orbital_space_b};

    // Zero-initialize a mixed spin tensor object.
    const auto M = orbital_space.initializeMixedRepresentableObjectFor<double>(occ, alpha, virt, beta, occ, beta, virt, alpha);

    // Test whether the dimensions of the mixed object match the correct number of indices of the chosen occupation type and spin.
    BOOST_CHECK(M.asTensor().dimension(0) == occupied_indices_a.size());
    BOOST_CHECK(M.asTensor().dimension(1) == virtual_indices_b.size());
    BOOST_CHECK(M.asTensor().dimension(2) == occupied_indices_b.size());
    BOOST_CHECK(M.asTensor().dimension(3) == virtual_indices_a.size());
}
