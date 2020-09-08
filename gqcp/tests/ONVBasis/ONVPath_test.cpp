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
BOOST_AUTO_TEST_CASE(sandbox) {

    // Set up a F(5,3) Fock space.
    const size_t M = 5;
    const size_t N = 3;

    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};


    // Set up the reference values from an example. [LINK TO WEBSITE]
    const auto I = GQCP::SpinUnresolvedONV::FromOccupiedIndices({0, 2, 3}, 5);  // |10110>
    GQCP::ONVPath onv_path {onv_basis, I};
    BOOST_REQUIRE(onv_path.address() == 2);

    // In the example, we're annihilating the first electron on the orbital with index 0. This means that we're annihilating the diagonal vertex that starts at (0,0).
    onv_path.annihilate(0, 0);

    // We've removed an arc with weight 0, so the current address hasn't changed.
    BOOST_REQUIRE(onv_path.address() == 2);


    // By manual inspection, we find that we may create a diagonal vertex on the current path, thereby closing it and forming a valid ONV.
    // After an annihilation, ONVPath sets its internal index that could be checked for a creation operator to (the previous annihilation index) + 1, which -in this case- is automatically correct.
    BOOST_REQUIRE(onv_path.nextCreationOrbitalIndex() == 1);


    // In order to close the current path, we should create a diagonal arc starting from the vertex (1,0).
    onv_path.create(1, 0);

    // We've managed to close the path, adding an arc with weight 1, so the address should now be 3.
    BOOST_REQUIRE(onv_path.address() == 3);


    // We'll have to undo the previous creation, in order to return to the path after the first annihilation.
    onv_path.annihilate(1, 0);

    // Since we've removed an arc with weight 1, the address should now be 2.
    BOOST_REQUIRE(onv_path.address() == 2);

    // The next orbital that should be checked for creation should now have index 2, since we previously annihilated on index 1).
    BOOST_REQUIRE(onv_path.nextCreationOrbitalIndex() == 2);


    // We're now in the situation that the next creation index corresponds to an occupied orbital. Since we can't create on these indices, we must translate the diagonal arc that starts at (2, 1) to (2, 0).
    onv_path.leftTranslate(2, 1);

    // A translation to the left means that we've encountered an electron, so the sign should be updated.
    BOOST_REQUIRE(onv_path.sign() == -1);

    // Since we've removed an arc with weight 1, and created one with weight 2, the address of the open path should now be 3.
    BOOST_REQUIRE(onv_path.address() == 3);


    // The next orbital that should be checked for creation should now have index 3, since we have moved up one orbital index.
    BOOST_REQUIRE(onv_path.nextCreationOrbitalIndex() == 3);


    // Since we're still in the situation that the next creation index corresponds to an occupied orbital, we repeat the previous procedure of left-translation.
    onv_path.leftTranslate(3, 2);

    // A translation to the left means that we've encountered an electron, so the sign should be updated.
    BOOST_REQUIRE(onv_path.sign() == 1);

    // Since we've removed an arc with weight 1, and created one with weight 3, the address of the open path should now be 5.
    BOOST_REQUIRE(onv_path.address() == 5);


    // The next orbital that should be checked for creation should now have index 3, since we have moved up one orbital index.
    BOOST_REQUIRE(onv_path.nextCreationOrbitalIndex() == 4);


    // We can now close up the path by creating an electron in the orbital with index 4 (the current creation index).
    onv_path.create(4, 2);

    // Since we've created an arc with weight 4, the address should now be 9.
    BOOST_REQUIRE(onv_path.address() == 9);

    // The total sign factor of this path is 1, since we've encountered 2 electrons in total.
    BOOST_REQUIRE(onv_path.sign() == 1);
}

/**
 *  Test if the shift in address, orbital indices and sign change correspond to the correct solutions
 *  
 *  This method does not alter the ONVs in any way.
 */
BOOST_AUTO_TEST_CASE(iterateToNextUnoccupiedOrbital) {

    // Set up a F(5,3) Fock space.
    const size_t M = 5;
    const size_t N = 3;

    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};


    // Set up the reference values from an example. [LINK TO WEBSITE]
    const auto I = GQCP::SpinUnresolvedONV::FromOccupiedIndices({0, 2, 3}, 5);  // |10110>
    GQCP::ONVPath onv_path {onv_basis, I};
    BOOST_REQUIRE(onv_path.address() == 2);

    // We start our example by annihilating the first electron on orbital index 0. This creates an open path.
    onv_path.annihilate(0, 0);

    // Starting from the index of the removed electron, we shift the "next possible creation operator" until we find an unoccupied orbital. In our case, this should be at p=1.
    onv_path.leftTranslateUntilVertical();
    const auto nextUnoccupiedIndex = onv_path.nextCreationOrbitalIndex();
    BOOST_REQUIRE(nextUnoccupiedIndex == 1);

    // We create the first electron up to this orbital.
    onv_path.create(nextUnoccupiedIndex, 0);

    // Subsequently, the second electron at orbital index 3 is annihilated.
    onv_path.annihilate(3, 1);

    // Starting from the removed electron, next unocupied orbital should be located at p=4.
    onv_path.leftTranslateUntilVertical();
    BOOST_REQUIRE(onv_path.nextCreationOrbitalIndex() == 4);





    //onv_path.shiftUntilNextUnoccupiedOrbital(2);
}
