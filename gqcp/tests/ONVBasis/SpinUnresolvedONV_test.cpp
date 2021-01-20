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

#define BOOST_TEST_MODULE "SpinUnresolvedONV"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "ONVBasis/SpinResolvedONV.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"


/**
 *  Check if the SpinUnresolvedONVBasis constructor is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(ONV_constructor) {

    // Create SpinUnresolvedONV. We consider in total 10 bits, with 5 set bits distributed as "0000011111" = 31
    GQCP::SpinUnresolvedONV onv1 {10, 5, 31};

    // Create SpinUnresolvedONV basis with 10 orbitals and 5 electrons, and check if its first ONV is the one we just constructed.
    const GQCP::SpinUnresolvedONVBasis onv_basis {10, 5};
    const auto onv2 = onv_basis.constructONVFromAddress(0);

    // Check if both unsigned representations match (are equal to 31 and have the same considered bits).
    BOOST_CHECK(onv1 == onv2);
    BOOST_CHECK(onv1.unsignedRepresentation() == 31);

    // Test if the occupation numbers are correct.
    std::vector<size_t> ref_indices {0, 1, 2, 3, 4};
    BOOST_CHECK(ref_indices == onv1.occupiedIndices());

    // Test if setting incompatible representation throws an error.
    BOOST_CHECK_THROW(onv1.replaceRepresentationWith(1), std::invalid_argument);
}


/**
 *  Check if the named constructor FromOccupiedIndices works as expected.
 */
BOOST_AUTO_TEST_CASE(FromOccupiedIndices) {

    const std::vector<size_t> occupied_indices1 {0, 2, 4};  // "10101" (21)
    const std::vector<size_t> occupied_indices2 {1, 2, 4};  // "10110" (22)
    const std::vector<size_t> occupied_indices3 {1, 3, 4};  // "11010" (26)

    const GQCP::SpinUnresolvedONV ref_onv1 {5, 3, 21};  // "10101" (21)
    const GQCP::SpinUnresolvedONV ref_onv2 {5, 3, 22};  // "10110" (22)
    const GQCP::SpinUnresolvedONV ref_onv3 {5, 3, 26};  // "11010" (26)

    const size_t M = 5;  // the total number of spinors
    BOOST_CHECK(GQCP::SpinUnresolvedONV::FromOccupiedIndices(occupied_indices1, M) == ref_onv1);
    BOOST_CHECK(GQCP::SpinUnresolvedONV::FromOccupiedIndices(occupied_indices2, M) == ref_onv2);
    BOOST_CHECK(GQCP::SpinUnresolvedONV::FromOccupiedIndices(occupied_indices3, M) == ref_onv3);
}


/**
 *  Check if the operator '==' behaves as expected.
 */
BOOST_AUTO_TEST_CASE(operator_equals_ONV) {

    const GQCP::SpinUnresolvedONV onv1 {6, 3, 19};  // "010011" (19)
    const GQCP::SpinUnresolvedONV onv2 {6, 3, 19};  // "010011" (19)
    const GQCP::SpinUnresolvedONV onv3 {7, 3, 19};  // "0010011" (19)

    BOOST_CHECK(onv1 == onv2);
    BOOST_CHECK(!(onv1 == onv3));  // The number of orbitals is wrong.
}


/**
 *  Check if the operator '!=' behaves as expected.
 */
BOOST_AUTO_TEST_CASE(operator_not_equals_ONV) {

    const GQCP::SpinUnresolvedONV onv1 {6, 3, 19};  // "010011" (19)
    const GQCP::SpinUnresolvedONV onv2 {6, 3, 19};  // "010011" (19)
    const GQCP::SpinUnresolvedONV onv3 {7, 3, 19};  // "0010011" (19)

    BOOST_CHECK(!(onv1 != onv2));
    BOOST_CHECK(onv1 != onv3);  // The number of orbitals is wrong.
}


/**
 *  Check if the method `isOccupied` behaves as expected.
 */
BOOST_AUTO_TEST_CASE(isOccupied) {

    const GQCP::SpinUnresolvedONV onv {4, 1, 2};  // "0010" (2)

    BOOST_CHECK_THROW(onv.isOccupied(9), std::invalid_argument);  // The index 9 is out of bounds.

    BOOST_CHECK(onv.isOccupied(1));
    BOOST_CHECK(!onv.isOccupied(2));
}


/**
 *  Check if the method `areOccupied` behaves as expected.
 */
BOOST_AUTO_TEST_CASE(areOccupied) {

    const GQCP::SpinUnresolvedONV onv {6, 3, 22};  // "010110" (22)

    BOOST_CHECK_THROW(onv.areOccupied({1, 2, 9}), std::invalid_argument);  // The index 9 is out of bounds.
    BOOST_CHECK_THROW(onv.areOccupied({9, 1, 2}), std::invalid_argument);  // The index 9 is out of bounds.

    BOOST_CHECK(onv.areOccupied({1, 2, 4}));
    BOOST_CHECK(onv.areOccupied({4, 2, 1}));

    BOOST_CHECK(!onv.areOccupied({0, 1}));
    BOOST_CHECK(!onv.areOccupied({3, 0}));
}


/**
 *  Check if the method `isUnoccupied` behaves as expected.
 */
BOOST_AUTO_TEST_CASE(isUnoccupied) {

    const GQCP::SpinUnresolvedONV onv {4, 1, 2};  // "0010" (2)

    BOOST_CHECK_THROW(onv.isUnoccupied(9), std::invalid_argument);  // The index 9 is out of bounds.

    BOOST_CHECK(onv.isUnoccupied(0));
    BOOST_CHECK(!onv.isUnoccupied(1));
}


/**
 *  Check if the method `areUnoccupied` behaves as expected.
 */
BOOST_AUTO_TEST_CASE(areUnoccupied) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};  // "010110" (22)

    BOOST_CHECK_THROW(onv.areUnoccupied({0, 3, 9}), std::invalid_argument);  // The index 9 is out of bounds.
    BOOST_CHECK_THROW(onv.areUnoccupied({9, 0, 3}), std::invalid_argument);  // The index 9 is out of bounds.

    BOOST_CHECK(onv.areUnoccupied({0, 3, 5}));
    BOOST_CHECK(onv.areUnoccupied({5, 0, 3}));

    BOOST_CHECK(!onv.areUnoccupied({0, 1}));
    BOOST_CHECK(!onv.areUnoccupied({1, 2}));
}


/**
 *  Check if the `slice` API behaves as expected.
 */
BOOST_AUTO_TEST_CASE(slice) {

    // Check if the error handling is correct.
    const GQCP::SpinUnresolvedONV onv1 {4, 2, 5};                // "0101" (5)
    BOOST_CHECK_THROW(onv1.slice(2, 1), std::invalid_argument);  // index_end should be larger than index_start
    BOOST_CHECK_THROW(onv1.slice(2, 2), std::invalid_argument);  // index_end should be larger than index_start
    BOOST_CHECK_THROW(onv1.slice(2, 6), std::invalid_argument);  // index_end is out of bounds


    // Check multiple test cases.
    const GQCP::SpinUnresolvedONV onv2 {7, 3, 41};
    BOOST_CHECK_EQUAL(onv2.slice(2, 6), 10);  // "0[1010]01" (41) -> "1010" (10)

    const GQCP::SpinUnresolvedONV onv3 {7, 5, 115};
    BOOST_CHECK_EQUAL(onv3.slice(2, 7), 28);  // "[11100]11" (115) -> "11100" (28)
    BOOST_CHECK_EQUAL(onv3.slice(2, 5), 4);   // "11[100]11" (115) -> "100" (4)

    const GQCP::SpinUnresolvedONV onv4 {6, 3, 19};
    BOOST_CHECK_EQUAL(onv4.slice(2, 6), 4);  // "[0100]11" (18) -> "0100" (4)
    BOOST_CHECK_EQUAL(onv4.slice(5, 6), 0);  // "[0]10011" (18) -> "0" (0)
    BOOST_CHECK_EQUAL(onv4.slice(4, 5), 1);  // "0[1]0011" (19) -> "1" (1)
}


/**
 *  Check if `operatorPhaseFactor` returns the correct sign: the sign should be negative on an index which has passed an odd number of electrons.
 */
BOOST_AUTO_TEST_CASE(operatorPhaseFactor) {

    // Test case 1.
    const GQCP::SpinUnresolvedONV onv1 {6, 3, 22};  // "010110" (22)
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(0), 1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(1), 1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(2), -1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(3), 1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(4), 1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(5), -1);

    // Test case 2.
    GQCP::SpinUnresolvedONV onv2 {6, 3, 26};  // "011010" (26)
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(0), 1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(1), 1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(2), -1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(3), -1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(4), 1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(5), -1);
}


/**
 *  Check if the method `annihilate` behaves as expected.
 */
BOOST_AUTO_TEST_CASE(annihilate) {

    GQCP::SpinUnresolvedONV onv {4, 2, 10};  // "1010" (10)

    // We shouldn't be able to annihilate on index 5 (out of bounds).
    BOOST_CHECK_THROW(onv.annihilate(5), std::invalid_argument);

    // We can't annihilate on index 2.
    BOOST_CHECK(!(onv.annihilate(2)));

    // We can annihilate on index 1.
    BOOST_CHECK(onv.annihilate(1));  // "1010" (10) -> "1000" (8)
    BOOST_CHECK_EQUAL(onv.unsignedRepresentation(), 8);

    // Test if updating throws an error (no longer 2 electrons).
    BOOST_CHECK_THROW(onv.updateOccupationIndices(), std::invalid_argument);
}


/**
 *  A first test to check the behavior of `annihilateAll`.
 */
BOOST_AUTO_TEST_CASE(annihilateAll_1) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};        // "010110" (22)
    const GQCP::SpinUnresolvedONV onv_copy = onv;  // Make a copy for further reference checks.

    BOOST_CHECK_THROW(onv.annihilateAll({8, 4}), std::invalid_argument);


    // Try to annihilate {0,1} and check if nothing changes.
    BOOST_CHECK(!onv.annihilateAll({0, 1}));
    BOOST_CHECK(onv_copy == onv);


    // Try to annihilate {1,0} and check if nothing changes.
    BOOST_CHECK(!onv.annihilateAll({1, 0}));
    BOOST_CHECK(onv_copy == onv);


    // Try to annihilate {5,2,0,1} and check if nothing changes.
    BOOST_CHECK(!onv.annihilateAll({5, 2, 0, 1}));
    BOOST_CHECK(onv_copy == onv);
}


/**
 *  A second test to check the behavior of `annihilateAll`.
 */
BOOST_AUTO_TEST_CASE(annihilateAll_2) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};  // "010110" (22)

    // Annihilate the indices {1, 2}.
    const GQCP::SpinUnresolvedONV ref_onv {6, 1, 16};  // "010000" (16)
    BOOST_CHECK(onv.annihilateAll({1, 2}));
    BOOST_CHECK(onv == ref_onv);
}


/**
 *  A third test to check the behavior of `annihilateAll`.
 */
BOOST_AUTO_TEST_CASE(annihilateAll_3) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};  // "010110" (22)

    // Annihilate the indices {2, 1}.
    const GQCP::SpinUnresolvedONV ref_onv {6, 1, 16};  // "010000" (16)
    BOOST_CHECK(onv.annihilateAll({2, 1}));
    BOOST_CHECK(onv == ref_onv);
}


/**
 *  Check if the `annihilate` method updates the sign correctly.
 */
BOOST_AUTO_TEST_CASE(annihilate_sign) {

    // There should be a sign change when we annihilate on (lexical) index 2 for "10101" (21).
    GQCP::SpinUnresolvedONV onv1 {5, 3, 21};
    int sign = 1;

    onv1.annihilate(2, sign);
    BOOST_CHECK_EQUAL(sign, -1);


    // There should be no sign change when we annihilate on (lexical) index 4 for "10101" (21).
    GQCP::SpinUnresolvedONV onv2 {5, 3, 21};
    sign = 1;

    onv2.annihilate(4, sign);
    BOOST_CHECK_EQUAL(sign, 1);


    // Annihilating the first occupied orbital on (lexical) index 0 for "10101" (21).
    GQCP::SpinUnresolvedONV onv3 {5, 3, 21};
    sign = 1;

    onv3.annihilate(0, sign);
    BOOST_CHECK_EQUAL(sign, 1);
}


/**
 *  A first test to check if the `annihilateAll` method updates the sign correctly.
 */
BOOST_AUTO_TEST_CASE(annihilateAll_sign_1) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};        // "010110" (22)
    const GQCP::SpinUnresolvedONV onv_copy = onv;  // Make a copy for further reference checks.
    int sign = 1;

    BOOST_CHECK_THROW(onv.annihilateAll({8, 4}, sign), std::invalid_argument);


    // Try to annihilate {0,1} and check if nothing changes.
    BOOST_CHECK(!onv.annihilateAll({0, 1}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);


    // Try to annihilate {1,0} and check if nothing changes.
    BOOST_CHECK(!onv.annihilateAll({1, 0}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);


    // Try to annihilate {5,2,0,1} and check if nothing changes.
    BOOST_CHECK(!onv.annihilateAll({5, 2, 0, 1}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);
}


/**
 *  A second test to check if the `annihilateAll` method updates the sign correctly.
 */
BOOST_AUTO_TEST_CASE(annihilateAll_sign_2) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};  // "010110" (22)

    // Annihilate the indices {1, 2}.
    const GQCP::SpinUnresolvedONV ref_onv {6, 1, 16};  // "010000" (16)
    int sign = 1;
    BOOST_CHECK(onv.annihilateAll({1, 2}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv == ref_onv);
}


/**
 *  A third test to check if the `annihilateAll` method updates the sign correctly.
 */
BOOST_AUTO_TEST_CASE(annihilateAll_sign_3) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};  // "010110" (22)

    // Annihilate the indices {2, 1}
    const GQCP::SpinUnresolvedONV ref_onv {6, 1, 16};  // "010000" (16)
    int sign = 1;
    BOOST_CHECK(onv.annihilateAll({2, 1}, sign));
    BOOST_CHECK(sign == -1);
    BOOST_CHECK(onv == ref_onv);
}


/**
 *  Check if the `create` method behaves as expected.
 */
BOOST_AUTO_TEST_CASE(create) {

    GQCP::SpinUnresolvedONV onv {4, 1, 2};  // "0010" (2)

    // We shouldn't be able to create on index 9 (out of bounds).
    BOOST_CHECK_THROW(onv.create(9), std::invalid_argument);

    // We can't create on index 1.
    BOOST_CHECK(!(onv.create(1)));

    // We can create on index 2.
    BOOST_CHECK(onv.create(2));  // "0010" (2) -> "0110" (6)
    BOOST_CHECK_EQUAL(onv.unsignedRepresentation(), 6);
}


/**
 *  A first test to check if the `createAll` method behaves correctly.
 */
BOOST_AUTO_TEST_CASE(createAll_1) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};        // "010110" (22)
    const GQCP::SpinUnresolvedONV onv_copy = onv;  // Make a copy for further reference checks.

    BOOST_CHECK_THROW(onv.createAll({8, 4}), std::invalid_argument);


    // Try to create {0,1} and check if nothing changes.
    BOOST_CHECK(!onv.createAll({0, 1}));
    BOOST_CHECK(onv_copy == onv);


    // Try to create {5,2,0,1} and check if nothing changes.
    BOOST_CHECK(!onv.createAll({5, 2, 0, 1}));
    BOOST_CHECK(onv_copy == onv);
}


/**
 *  A second test to check if the `createAll` method updates the sign correctly.
 */
BOOST_AUTO_TEST_CASE(createAll_2) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};  // "010110" (22)

    // Create the indices {0, 3}.
    const GQCP::SpinUnresolvedONV ref_onv {6, 5, 31};  // "011111" (31)
    BOOST_CHECK(onv.createAll({0, 3}));
    BOOST_CHECK(onv == ref_onv);
}


/**
 *  Check if the `create` method updates the sign correctly.
 */
BOOST_AUTO_TEST_CASE(create_sign) {

    // There should be a sign change when we create on (lexical) index 1 for "10101" (21).
    GQCP::SpinUnresolvedONV onv1 {5, 3, 21};
    int sign = 1;

    onv1.create(1, sign);
    BOOST_CHECK_EQUAL(sign, -1);


    // There should be no sign change when we create on (lexical) index 3 for "10101" (21).
    GQCP::SpinUnresolvedONV onv2 {5, 3, 21};
    sign = 1;

    onv2.create(3, sign);
    BOOST_CHECK_EQUAL(sign, 1);


    // Creating the last orbital on (lexical) index 4 for "00101" (5).
    GQCP::SpinUnresolvedONV onv3 {5, 2, 5};
    sign = 1;

    onv3.create(4, sign);
    BOOST_CHECK_EQUAL(sign, 1);


    // Creating the first orbital on (lexical) index 0 for "10100" (5).
    GQCP::SpinUnresolvedONV onv4 {5, 2, 20};
    sign = 1;

    onv4.create(0, sign);
    BOOST_CHECK_EQUAL(sign, 1);
}


/**
 *  A first test to check if the `createAll` method updates the sign correctly.
 */
BOOST_AUTO_TEST_CASE(createAll_sign_1) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};        // "010110" (22)
    const GQCP::SpinUnresolvedONV onv_copy = onv;  // Make a copy for further reference checks.

    BOOST_CHECK_THROW(onv.createAll({8, 4}), std::invalid_argument);


    // Try to create {0,1} and check if nothing changes.
    int sign = 1;
    BOOST_CHECK(!onv.createAll({0, 1}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);


    // Try to create {5,2,0,1} and check if nothing changes.
    BOOST_CHECK(!onv.createAll({5, 2, 0, 1}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);
}


/**
 *  A second test to check if the `createAll` method updates the sign correctly.
 */
BOOST_AUTO_TEST_CASE(createAll_sign_2) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};  // "010110" (22)

    // Create the indices {0, 3}.
    const GQCP::SpinUnresolvedONV ref_onv {6, 5, 31};  // "011111" (31)
    int sign = 1;
    BOOST_CHECK(onv.createAll({0, 3}, sign));
    BOOST_CHECK(sign == -1);
    BOOST_CHECK(onv == ref_onv);
}


/**
 *  A third test to check if the `createAll` method updates the sign correctly.
 */
BOOST_AUTO_TEST_CASE(createAll_sign_3) {

    GQCP::SpinUnresolvedONV onv {6, 3, 22};  // "010110" (22)

    // Create the indices {3, 0}
    int sign = 1;
    const GQCP::SpinUnresolvedONV ref_onv {6, 5, 31};  // "011111" (31)
    BOOST_CHECK(onv.createAll({3, 0}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv == ref_onv);
}


/**
 *  Check if the `countNumberOfDifferences` method behaves correctly.
 */
BOOST_AUTO_TEST_CASE(countNumberOfDifferences) {

    const GQCP::SpinUnresolvedONV onv1 {5, 3, 21};  // "10101" (21)
    const GQCP::SpinUnresolvedONV onv2 {5, 3, 22};  // "10110" (22)
    const GQCP::SpinUnresolvedONV onv3 {5, 3, 26};  // "11010" (26)

    BOOST_CHECK_EQUAL(onv1.countNumberOfDifferences(onv1), 0);
    BOOST_CHECK_EQUAL(onv2.countNumberOfDifferences(onv2), 0);
    BOOST_CHECK_EQUAL(onv3.countNumberOfDifferences(onv3), 0);

    BOOST_CHECK_EQUAL(onv1.countNumberOfDifferences(onv2), 2);
    BOOST_CHECK_EQUAL(onv1.countNumberOfDifferences(onv3), 4);
    BOOST_CHECK_EQUAL(onv2.countNumberOfDifferences(onv3), 2);
}


/**
 *  Check if the `countNumberOfDifferences` method behaves correctly.
 */
BOOST_AUTO_TEST_CASE(findDifferentOccupations) {

    const GQCP::SpinUnresolvedONV onv1 {5, 3, 21};  // "10101" (21)
    const GQCP::SpinUnresolvedONV onv2 {5, 3, 22};  // "10110" (22)
    const GQCP::SpinUnresolvedONV onv3 {5, 3, 26};  // "11010" (26)

    BOOST_TEST(onv1.findDifferentOccupations(onv2) == (std::vector<size_t> {0}), boost::test_tools::per_element());
    BOOST_TEST(onv2.findDifferentOccupations(onv1) == (std::vector<size_t> {1}), boost::test_tools::per_element());

    BOOST_TEST(onv1.findDifferentOccupations(onv3) == (std::vector<size_t> {0, 2}), boost::test_tools::per_element());
    BOOST_TEST(onv3.findDifferentOccupations(onv1) == (std::vector<size_t> {1, 3}), boost::test_tools::per_element());

    BOOST_TEST(onv2.findDifferentOccupations(onv3) == (std::vector<size_t> {2}), boost::test_tools::per_element());
    BOOST_TEST(onv3.findDifferentOccupations(onv2) == (std::vector<size_t> {3}), boost::test_tools::per_element());
}


/**
 *  Check if the `findMatchingOccupations` method behaves correctly.
 */
BOOST_AUTO_TEST_CASE(findMatchingOccupations) {

    const GQCP::SpinUnresolvedONV onv1 {5, 3, 21};  // "10101" (21)
    const GQCP::SpinUnresolvedONV onv2 {5, 3, 22};  // "10110" (22)
    const GQCP::SpinUnresolvedONV onv3 {5, 3, 26};  // "11010" (26)

    BOOST_TEST(onv1.findMatchingOccupations(onv2) == (std::vector<size_t> {2, 4}), boost::test_tools::per_element());
    BOOST_TEST(onv2.findMatchingOccupations(onv1) == (std::vector<size_t> {2, 4}), boost::test_tools::per_element());

    BOOST_TEST(onv1.findMatchingOccupations(onv3) == (std::vector<size_t> {4}), boost::test_tools::per_element());
    BOOST_TEST(onv3.findMatchingOccupations(onv1) == (std::vector<size_t> {4}), boost::test_tools::per_element());

    BOOST_TEST(onv2.findMatchingOccupations(onv3) == (std::vector<size_t> {1, 4}), boost::test_tools::per_element());
    BOOST_TEST(onv3.findMatchingOccupations(onv2) == (std::vector<size_t> {1, 4}), boost::test_tools::per_element());
}


/**
 *  Check if the creation of the 'GHF' ONV works as expected.
 */
BOOST_AUTO_TEST_CASE(GHF) {

    // Create an example array of orbital energies.
    const size_t M = 6;  // The total number of spinors.
    GQCP::VectorX<double> orbital_energies {M};
    orbital_energies << 0, 1, 2, 0, 0.5, 1;

    // N = 2.
    const auto ref_onv1 = GQCP::SpinUnresolvedONV::FromString("001001");
    BOOST_CHECK(GQCP::SpinUnresolvedONV::GHF(M, 2, orbital_energies) == ref_onv1);

    // N = 3.
    const auto ref_onv2 = GQCP::SpinUnresolvedONV::FromString("011001");
    BOOST_CHECK(GQCP::SpinUnresolvedONV::GHF(M, 3, orbital_energies) == ref_onv2);

    // N = 4.
    const auto ref_onv3 = GQCP::SpinUnresolvedONV::FromString("011011");
    BOOST_CHECK(GQCP::SpinUnresolvedONV::GHF(M, 4, orbital_energies) == ref_onv3);

    // N = 5.
    const auto ref_onv4 = GQCP::SpinUnresolvedONV::FromString("111011");
    BOOST_CHECK(GQCP::SpinUnresolvedONV::GHF(M, 5, orbital_energies) == ref_onv4);
}


/**
 *  Check if the creation of a spin-unresolved ONV from a textual/string representation works as expected.
 */
BOOST_AUTO_TEST_CASE(FromString) {

    const GQCP::SpinUnresolvedONV onv_ref1 {5, 3, 21};  // "10101" (21)
    const GQCP::SpinUnresolvedONV onv_ref2 {5, 3, 22};  // "10110" (22)
    const GQCP::SpinUnresolvedONV onv_ref3 {5, 3, 26};  // "11010" (26)

    const auto onv1 = GQCP::SpinUnresolvedONV::FromString("10101");
    BOOST_CHECK(onv1 == onv_ref1);

    const auto onv2 = GQCP::SpinUnresolvedONV::FromString("10110");
    BOOST_CHECK(onv2 == onv_ref2);

    const auto onv3 = GQCP::SpinUnresolvedONV::FromString("11010");
    BOOST_CHECK(onv3 == onv_ref3);
}


/**
 *  Check if unoccupiedIndices() works as expected.
 */
BOOST_AUTO_TEST_CASE(unoccupiedIndices) {

    // Create some test ONVs.
    const GQCP::SpinUnresolvedONV onv1 {5, 3, 21};  // "10101" (21)
    const GQCP::SpinUnresolvedONV onv2 {5, 3, 22};  // "10110" (22)
    const GQCP::SpinUnresolvedONV onv3 {5, 3, 26};  // "11010" (26)
    const GQCP::SpinUnresolvedONV onv4 {4, 2, 5};   // "0101" (5)

    // Initialize the reference unoccupied indices and check the result.
    const std::vector<size_t> ref_unoccupied_indices1 {1, 3};
    const std::vector<size_t> ref_unoccupied_indices2 {0, 3};
    const std::vector<size_t> ref_unoccupied_indices3 {0, 2};
    const std::vector<size_t> ref_unoccupied_indices4 {1, 3};

    BOOST_CHECK(onv1.unoccupiedIndices() == ref_unoccupied_indices1);  // Counted from right to left.
    BOOST_CHECK(onv2.unoccupiedIndices() == ref_unoccupied_indices2);  // Counted from right to left.
    BOOST_CHECK(onv3.unoccupiedIndices() == ref_unoccupied_indices3);  // Counted from right to left.
    BOOST_CHECK(onv4.unoccupiedIndices() == ref_unoccupied_indices4);  // Counted from right to left.
}


/**
 *  Check if the overlap between a two GHF-related ONVs works as expected.
 *
 *  We don't really have a reference implementation, but we can check if the overlaps are equal to 1 or 0 if we use GHF orbitals (constructed from RHF canonical spin-orbitals) that are equal.
 *
 *  The system under consideration is H2 with a STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(GHF_overlap) {

    // Obtain the canonical RHF spin-orbitals.
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const auto N = h2.numberOfElectrons();
    const auto N_P = h2.numberOfElectronPairs();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> r_spinor_basis {h2, "STO-3G"};
    const auto K = r_spinor_basis.numberOfSpatialOrbitals();
    const auto S_restricted = r_spinor_basis.overlap().parameters();

    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(r_spinor_basis, h2);  // In an AO basis.

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, S_restricted);
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};

    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    // Create a GSpinorBasis from the canonical RHF spin-orbitals, yielding a general spinor basis with alternating alpha- and beta-spin-orbitals.
    r_spinor_basis.transform(rhf_parameters.expansion());

    const auto g_spinor_basis_of = GQCP::GSpinorBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);
    const auto g_spinor_basis_on = g_spinor_basis_of;

    const auto& C_on = g_spinor_basis_on.expansion();
    const auto& C_of = g_spinor_basis_of.expansion();

    auto S_op = g_spinor_basis_of.overlap();  // In the GHF MO basis.
    S_op.transform(C_on.inverse());           // Now in back in the AO basis.
    const auto& S = S_op.parameters();


    // Check if the one GHF determinant has overlap 1 with the other corresponding GHF determinant, and overlap 0 with the other excitations.
    // For the construction of the GHF determinant, we need the spinor orbital energies.
    GQCP::VectorX<double> total_orbital_energies {2 * K};
    total_orbital_energies.topRows(K) = rhf_parameters.orbitalEnergies();
    total_orbital_energies.bottomRows(K) = rhf_parameters.orbitalEnergies();

    const auto ghf_determinant_of = GQCP::SpinUnresolvedONV::GHF(2 * K, N, total_orbital_energies);


    // Create all the excitations and check if the projection gives the right value.
    const auto onv_0011 = GQCP::SpinUnresolvedONV::FromString("0011");
    const auto onv_0101 = GQCP::SpinUnresolvedONV::FromString("0101");  // this is the GHF determinant (according to the orbital energies)
    const auto onv_0110 = GQCP::SpinUnresolvedONV::FromString("0110");
    const auto onv_1001 = GQCP::SpinUnresolvedONV::FromString("1001");
    const auto onv_1010 = GQCP::SpinUnresolvedONV::FromString("1010");
    const auto onv_1100 = GQCP::SpinUnresolvedONV::FromString("1100");

    BOOST_CHECK(std::abs(ghf_determinant_of.calculateProjection(onv_0011, C_on, C_of, S) - 0.0) < 1.0e-12);
    BOOST_CHECK(std::abs(ghf_determinant_of.calculateProjection(onv_0101, C_on, C_of, S) - 1.0) < 1.0e-12);
    BOOST_CHECK(std::abs(ghf_determinant_of.calculateProjection(onv_0110, C_on, C_of, S) - 0.0) < 1.0e-12);
    BOOST_CHECK(std::abs(ghf_determinant_of.calculateProjection(onv_1001, C_on, C_of, S) - 0.0) < 1.0e-12);
    BOOST_CHECK(std::abs(ghf_determinant_of.calculateProjection(onv_1010, C_on, C_of, S) - 0.0) < 1.0e-12);
    BOOST_CHECK(std::abs(ghf_determinant_of.calculateProjection(onv_1100, C_on, C_of, S) - 0.0) < 1.0e-12);
}


/**
 *  Check if projecting <RHF|UHF> can also be calculated from <GHF|GHF>, i.e. check if the specialized algorithm is consistent with the general algorithm.
 * 
 *  The system of interest is an H4-square in an STO-3G basisset. The RHF and UHF results were found using Xeno's GHF code.
 */
BOOST_AUTO_TEST_CASE(RHF_UHF_projection) {

    const auto molecule = GQCP::Molecule::HRingFromDistance(4, 1.0);

    // Obtain the canonical RHF spin-orbitals.
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
    const auto S_restricted = r_spinor_basis.overlap().parameters();  // The AO overlap matrix.
    GQCP::SquareMatrix<double> C_restricted_matrix {4};               // The RHF canonical orbitals for this system, from a prototype Python calculation.
    // clang-format off
    C_restricted_matrix << -0.27745359, -0.8505133,   0.85051937,  2.02075317,
                           -0.27745362, -0.85051937, -0.8505133,  -2.02075317,
                           -0.27745359,  0.8505133,  -0.85051937,  2.02075317,
                           -0.27745362,  0.85051937,  0.8505133,  -2.02075317;
    // clang-format on
    GQCP::RTransformation<double> C_restricted {C_restricted_matrix};
    r_spinor_basis.transform(C_restricted);


    // Obtain the canonical UHF spin-orbitals.
    GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> u_spinor_basis {molecule, "STO-3G"};
    GQCP::SquareMatrix<double> C_alpha_matrix {4};  // The UHF alpha canonical orbitals for this system, from a prototype Python calculation.
    // clang-format off
    C_alpha_matrix << -1.75646828e-01, -1.20606646e-06,  1.20281173e+00,  2.03213486e+00,
                      -3.78560533e-01, -1.20281173e+00, -1.20606647e-06, -2.00427438e+00,
                      -1.75646828e-01,  1.20606646e-06, -1.20281173e+00,  2.03213486e+00,
                      -3.78560533e-01,  1.20281173e+00,  1.20606646e-06, -2.00427438e+00;
    // clang-format on

    GQCP::SquareMatrix<double> C_beta_matrix {4};  // The UHF beta canonical orbitals for this system, from a prototype Python calculation.
    // clang-format off
    C_beta_matrix << -3.78560533e-01,  1.20281173e+00,  1.21724557e-06,  2.00427438e+00,
                     -1.75646828e-01,  1.21724558e-06, -1.20281173e+00, -2.03213486e+00,
                     -3.78560533e-01, -1.20281173e+00, -1.21724558e-06,  2.00427438e+00,
                     -1.75646828e-01, -1.21724558e-06,  1.20281173e+00, -2.03213486e+00;
    // clang-format on
    GQCP::UTransformation<double> C_unrestricted {C_alpha_matrix, C_beta_matrix};
    u_spinor_basis.transform(C_unrestricted);


    // Calculate <RHF|UHF> using the specialized formula.
    const auto rhf_determinant = GQCP::SpinResolvedONV::RHF(4, 2);     // 4 spatial orbitals, 2 electron pairs.
    const auto uhf_determinant = GQCP::SpinResolvedONV::UHF(4, 2, 2);  // 4 spatial orbitals for each spin component, 2 alpha electrons, 2 beta electrons.

    const auto uhf_on_rhf_projection_specialized = uhf_determinant.calculateProjection(rhf_determinant, C_unrestricted, C_restricted, S_restricted);


    // Calculate <RHF|UHF> using the general formula.
    const auto g_spinor_basis_of = GQCP::GSpinorBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);
    const auto g_spinor_basis_on = GQCP::GSpinorBasis<double, GQCP::GTOShell>::FromUnrestricted(u_spinor_basis);

    const auto& C_of = g_spinor_basis_of.expansion();
    const auto& C_on = g_spinor_basis_on.expansion();

    auto S_generalized_op = g_spinor_basis_of.overlap();  // In MO basis.
    S_generalized_op.transform(C_of.inverse());           // In AO basis.
    const auto& S_generalized = S_generalized_op.parameters();


    const auto onv_on = GQCP::SpinUnresolvedONV::FromString("00110011");
    const auto onv_of = GQCP::SpinUnresolvedONV::FromString("00110011");

    const auto uhf_on_rhf_projection_general = onv_of.calculateProjection(onv_on, C_of, C_on, S_generalized);


    // Check if both approaches yield the same result.
    BOOST_CHECK(std::abs(uhf_on_rhf_projection_general - uhf_on_rhf_projection_specialized) < 1.0e-12);
}
