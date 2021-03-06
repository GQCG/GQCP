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

#define BOOST_TEST_MODULE "ONVPath"

#include <boost/test/unit_test.hpp>

#include "ONVBasis/ONVPath.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"


/**
 *  Check if we can use ONVPath to find ONVs that are one electron excitation away from a reference.
 */
BOOST_AUTO_TEST_CASE(example) {

    // Set up a F(5,3) Fock space.
    const size_t M = 5;
    const size_t N = 3;

    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};


    // Set up the reference values from an example. [https://gqcg.github.io/GQCP/docs/developer_documentation/ONV_path_manipulation]
    const auto I = GQCP::SpinUnresolvedONV::FromOccupiedIndices({0, 2, 3}, 5);  // |10110>
    GQCP::ONVPath<GQCP::SpinUnresolvedONVBasis> onv_path {onv_basis, I};
    BOOST_REQUIRE(onv_path.address() == 2);

    // In the example, we're annihilating the first electron on the orbital with index 0. This means that we're annihilating the diagonal vertex that starts at (0,0).
    onv_path.annihilate(0, 0);

    // We've removed an arc with weight 0, so the current address hasn't changed.
    BOOST_REQUIRE(onv_path.address() == 2);


    // By manual inspection, we find that we may create a diagonal vertex on the current path, thereby closing it and forming a valid ONV.
    // After an annihilation, ONVPath sets its internal index that could be checked for a creation operator to (the previous annihilation index) + 1, which -in this case- is automatically correct.
    BOOST_REQUIRE(onv_path.orbitalIndex() == 1);


    // In order to close the current path, we should create a diagonal arc starting from the vertex (1,0).
    onv_path.create(1, 0);

    // We've managed to close the path, adding an arc with weight 1, so the address should now be 3.
    BOOST_REQUIRE(onv_path.address() == 3);


    // We'll have to undo the previous creation, in order to return to the path after the first annihilation.
    onv_path.annihilate(1, 0);

    // Since we've removed an arc with weight 1, the address should now be 2.
    BOOST_REQUIRE(onv_path.address() == 2);

    // The next orbital that should be checked for creation should now have index 2, since we previously annihilated on index 1).
    BOOST_REQUIRE(onv_path.orbitalIndex() == 2);


    // We're now in the situation that the next creation index corresponds to an occupied orbital. Since we can't create on these indices, we must translate the diagonal arc that starts at (2, 1) to (2, 0).
    onv_path.leftTranslateDiagonalArc(2, 1);

    // A translation to the left means that we've encountered an electron, so the sign should be updated.
    BOOST_REQUIRE(onv_path.sign() == -1);

    // Since we've removed an arc with weight 1, and created one with weight 2, the address of the open path should now be 3.
    BOOST_REQUIRE(onv_path.address() == 3);


    // The next orbital that should be checked for creation should now have index 3, since we have moved up one orbital index.
    BOOST_REQUIRE(onv_path.orbitalIndex() == 3);


    // Since we're still in the situation that the next creation index corresponds to an occupied orbital, we repeat the previous procedure of left-translation.
    onv_path.leftTranslateDiagonalArc(3, 2);

    // A translation to the left means that we've encountered an electron, so the sign should be updated.
    BOOST_REQUIRE(onv_path.sign() == 1);

    // Since we've removed an arc with weight 1, and created one with weight 3, the address of the open path should now be 5.
    BOOST_REQUIRE(onv_path.address() == 5);


    // The next orbital that should be checked for creation should now have index 3, since we have moved up one orbital index.
    BOOST_REQUIRE(onv_path.orbitalIndex() == 4);


    // We can now close up the path by creating an electron in the orbital with index 4 (the current creation index).
    onv_path.create(4, 2);

    // Since we've created an arc with weight 4, the address should now be 9.
    BOOST_REQUIRE(onv_path.address() == 9);

    // The total sign factor of this path is 1, since we've encountered 2 electrons in total.
    BOOST_REQUIRE(onv_path.sign() == 1);
}


/**
 *  Check if the annihilate function of the ONV path works as expected.
 */
BOOST_AUTO_TEST_CASE(annihilate) {

    // Set up a F(5,3) Fock space.
    const size_t M = 5;
    const size_t N = 3;

    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};


    // Set up the reference values from an example. [https://gqcg.github.io/GQCP/docs/developer_documentation/ONV_path_manipulation]
    const auto I = GQCP::SpinUnresolvedONV::FromOccupiedIndices({0, 2, 3}, 5);  // |10110>
    GQCP::ONVPath<GQCP::SpinUnresolvedONVBasis> onv_path {onv_basis, I};
    GQCP::ONVPath<GQCP::SpinUnresolvedONVBasis> onv_path_2 {onv_basis, I};
    BOOST_REQUIRE(onv_path.address() == onv_path_2.address() && onv_path.address() == 2);

    // Current path state should at this point be the same as (0, 0)
    onv_path.annihilate();
    onv_path_2.annihilate(0, 0);

    BOOST_REQUIRE(onv_path.address() == onv_path_2.address() && onv_path.address() == 2);

    onv_path.create();
    onv_path_2.create();

    onv_path.annihilate();
    onv_path_2.annihilate(1, 0);

    BOOST_REQUIRE(onv_path.address() == onv_path_2.address() && onv_path.address() == 2);
}


/**
 *  Check if the create function of the ONV path works as expected.
 */
BOOST_AUTO_TEST_CASE(create) {

    // Set up a F(5,3) Fock space.
    const size_t M = 5;
    const size_t N = 3;

    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};


    // Set up the reference values from an example. [https://gqcg.github.io/GQCP/docs/developer_documentation/ONV_path_manipulation]
    const auto I = GQCP::SpinUnresolvedONV::FromOccupiedIndices({0, 2, 3}, 5);  // |10110>
    GQCP::ONVPath<GQCP::SpinUnresolvedONVBasis> onv_path {onv_basis, I};
    GQCP::ONVPath<GQCP::SpinUnresolvedONVBasis> onv_path_2 {onv_basis, I};
    BOOST_REQUIRE(onv_path.address() == onv_path_2.address() && onv_path.address() == 2);

    onv_path.annihilate();
    onv_path_2.annihilate();

    // Current path state should at this point be the same as (1, 0)
    onv_path.create();
    onv_path_2.create(1, 0);

    BOOST_REQUIRE(onv_path.address() == onv_path_2.address() && onv_path.address() == 3);
}


/**
 *  Test if the shift in address, orbital indices and sign change correspond to the correct solutions.
 */
BOOST_AUTO_TEST_CASE(leftTranslateUntilVerticalArc) {

    // Set up a F(5,3) Fock space.
    const size_t M = 5;
    const size_t N = 3;

    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};


    // Set up the reference values from an example. [https://gqcg.github.io/GQCP/docs/developer_documentation/ONV_path_manipulation]
    const auto I = GQCP::SpinUnresolvedONV::FromOccupiedIndices({1, 2, 3}, 5);  // |01110>
    GQCP::ONVPath<GQCP::SpinUnresolvedONVBasis> onv_path {onv_basis, I};
    GQCP::ONVPath<GQCP::SpinUnresolvedONVBasis> onv_path_ref {onv_basis, I};
    BOOST_REQUIRE(onv_path.address() == onv_path_ref.address() && onv_path.address() == 3);

    // Since we have removed a diagonal at position (1, 0) with arc weight 1, the address will change.
    onv_path.annihilate(1, 0);
    onv_path_ref.annihilate(1, 0);
    BOOST_REQUIRE(onv_path.address() == onv_path_ref.address() && onv_path.address() == 2);

    // To find the next creation index, i.e. a vertical arc sign signifying an unoccupied orbital, we must translate the diagonal arcs one to the left until we find a vertical arc.
    onv_path_ref.leftTranslateDiagonalArc(2, 1);

    // The diagonal arc at (2, 1) with arc weight = 1 is moved to (2,0) with arc weight = 2 so the address will change.
    BOOST_REQUIRE(onv_path_ref.address() == 3);

    // The next orbital index at (3, 2) is occupied, so we must translate it likewise.
    onv_path_ref.leftTranslateDiagonalArc(3, 2);

    // The arc weight at (3, 2) is 1 while the arc weight at the new occupation index (3, 1) is 3, meaning the address will change.
    BOOST_REQUIRE(onv_path_ref.address() == 5);

    // We have now found our next creation index, located at (4,2).
    BOOST_REQUIRE(onv_path_ref.orbitalIndex() == 4 && onv_path_ref.electronIndex() == 2);

    // This manual result must be the same as the result obtained by leftTranslateUntilVertical().
    onv_path.leftTranslateDiagonalArcUntilVerticalArc();
    BOOST_REQUIRE(onv_path.orbitalIndex() == onv_path_ref.orbitalIndex() && onv_path.electronIndex() == onv_path_ref.electronIndex());
}
