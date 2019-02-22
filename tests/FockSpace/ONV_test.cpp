// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#define BOOST_TEST_MODULE "ONV"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain

#include "FockSpace/ONV.hpp"
#include "FockSpace/FockSpace.hpp"



BOOST_AUTO_TEST_CASE ( ONV_constructor ) {

    // Create ONV : 10 considered bits and 5 set bits with distributed as "0000011111" = 31
    GQCP::ONV onv1 (10, 5, 31);

    // Create Fock space with 10 orbitals and 5 electrons
    GQCP::FockSpace fock_space (10, 5);

    // Ask for the first ONV in reverse lexicographic order : "0000011111" = 31
    GQCP::ONV onv2 = fock_space.makeONV(0);

    // Check if both unsigned representations match (are equal to 31 and have the same considered bits)
    BOOST_CHECK(onv1 == onv2);
    BOOST_CHECK(onv1.get_unsigned_representation() == 31);

    // Test if the occupation numbers are correct.
    GQCP::VectorXs x (5);
    x << 0, 1, 2, 3, 4;
    BOOST_CHECK(x.isApprox(onv1.get_occupation_indices()));

    // Test if setting incompatible representation throws an error
    BOOST_CHECK_THROW(onv1.set_representation(1), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( operator_equals_ONV ) {

    GQCP::ONV onv1 (6, 3, 19);  // "010011" (19)
    GQCP::ONV onv2 (6, 3, 19);  // "010011" (19)
    GQCP::ONV onv3 (7, 3, 19);  // "0010011" (19)

    BOOST_CHECK(onv1 == onv2);
    BOOST_CHECK(!(onv1 == onv3));  // wrong number of orbitals
}


BOOST_AUTO_TEST_CASE ( operator_not_equals_ONV ) {

    GQCP::ONV onv1 (6, 3, 19);  // "010011" (19)
    GQCP::ONV onv2 (6, 3, 19);  // "010011" (19)
    GQCP::ONV onv3 (7, 3, 19);  // "0010011" (19)

    BOOST_CHECK(!(onv1 != onv2));
    BOOST_CHECK(onv1 != onv3);  // wrong number of orbitals
}


BOOST_AUTO_TEST_CASE ( isOccupied ) {

    GQCP::ONV onv (4, 1, 2);  // "0010" (2)

    BOOST_CHECK_THROW(onv.isOccupied(9), std::invalid_argument);  // 9 is out of bounds

    BOOST_CHECK(onv.isOccupied(1));
    BOOST_CHECK(!onv.isOccupied(2));
}


BOOST_AUTO_TEST_CASE ( areOccupied ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)

    BOOST_CHECK_THROW(onv.areOccupied({1, 2, 9}), std::invalid_argument);  // 9 is out of bounds
    BOOST_CHECK_THROW(onv.areOccupied({9, 1, 2}), std::invalid_argument);  // 9 is out of bounds

    BOOST_CHECK(onv.areOccupied({1, 2, 4}));
    BOOST_CHECK(onv.areOccupied({4, 2, 1}));

    BOOST_CHECK(!onv.areOccupied({0, 1}));
    BOOST_CHECK(!onv.areOccupied({3, 0}));
}


BOOST_AUTO_TEST_CASE ( isUnoccupied ) {

    GQCP::ONV onv (4, 1, 2);  // "0010" (2)

    BOOST_CHECK_THROW(onv.isUnoccupied(9), std::invalid_argument);  // 9 is out of bounds

    BOOST_CHECK(onv.isUnoccupied(0));
    BOOST_CHECK(!onv.isUnoccupied(1));
}


BOOST_AUTO_TEST_CASE ( areUnoccupied ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)

    BOOST_CHECK_THROW(onv.areUnoccupied({0, 3, 9}), std::invalid_argument);  // 9 is out of bounds
    BOOST_CHECK_THROW(onv.areUnoccupied({9, 0, 3}), std::invalid_argument);  // 9 is out of bounds

    BOOST_CHECK(onv.areUnoccupied({0, 3, 5}));
    BOOST_CHECK(onv.areUnoccupied({5, 0, 3}));

    BOOST_CHECK(!onv.areUnoccupied({0, 1}));
    BOOST_CHECK(!onv.areUnoccupied({1, 2}));
}


BOOST_AUTO_TEST_CASE ( slice ) {

    // Check for throws
    GQCP::ONV onv1 (4, 2, 5);  // "0101" (5)
    BOOST_CHECK_THROW(onv1.slice(2, 1), std::invalid_argument);  // index_end should be larger than index_start
    BOOST_CHECK_THROW(onv1.slice(2, 2), std::invalid_argument);  // index_end should be larger than index_start
    BOOST_CHECK_THROW(onv1.slice(2, 6), std::invalid_argument);  // index_end is out of bounds


    GQCP::ONV onv2 (7, 3, 41);
    BOOST_CHECK_EQUAL(onv2.slice(2, 6), 10);  // "0[1010]01" (41) -> "1010" (10)

    GQCP::ONV onv3 (7, 5, 115);
    BOOST_CHECK_EQUAL(onv3.slice(2, 7), 28);  // "[11100]11" (115) -> "11100" (28)
    BOOST_CHECK_EQUAL(onv3.slice(2, 5), 4);   // "11[100]11" (115) -> "100" (4)

    GQCP::ONV onv4 (6, 3, 19);
    BOOST_CHECK_EQUAL(onv4.slice(2, 6), 4);  // "[0100]11" (18) -> "0100" (4)
    BOOST_CHECK_EQUAL(onv4.slice(5, 6), 0);  // "[0]10011" (18) -> "0" (0)
    BOOST_CHECK_EQUAL(onv4.slice(4, 5), 1);  // "0[1]0011" (19) -> "1" (1)
}


BOOST_AUTO_TEST_CASE ( operatorPhaseFactor ) {

    // The sign should be negative on an index which has passed an odd amount of electrons
    GQCP::ONV onv1 (6, 3, 22);  // "010110" (22)

    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(0), 1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(1), 1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(2), -1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(3), 1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(4), 1);
    BOOST_CHECK_EQUAL(onv1.operatorPhaseFactor(5), -1);

    GQCP::ONV onv2 (6, 3, 26);  // "011010" (26)
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(0), 1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(1), 1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(2), -1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(3), -1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(4), 1);
    BOOST_CHECK_EQUAL(onv2.operatorPhaseFactor(5), -1);
}


BOOST_AUTO_TEST_CASE ( annihilate ) {

    GQCP::ONV onv (4, 2, 10);  // "1010" (10)

    // We shouldn't be able to annihilate on index 5 (out of bounds)
    BOOST_CHECK_THROW(onv.annihilate(5), std::invalid_argument);

    // We can't annihilate on index 2
    BOOST_CHECK(!(onv.annihilate(2)));

    // We can annihilate on index 1
    BOOST_CHECK(onv.annihilate(1));  // "1010" (10) -> "1000" (8)
    BOOST_CHECK_EQUAL(onv.get_unsigned_representation(), 8);
    // Test if updating throws an error (no longer 2 electrons)
    BOOST_CHECK_THROW(onv.updateOccupationIndices(), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( annihilateAll_1 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)
    GQCP::ONV onv_copy = onv;  // to check if nothing changes

    BOOST_CHECK_THROW(onv.annihilateAll({8, 4}), std::invalid_argument);


    // Try to annihilate {0,1} and check if nothing changes
    BOOST_CHECK(!onv.annihilateAll({0, 1}));
    BOOST_CHECK(onv_copy == onv);


    // Try to annihilate {1,0} and check if nothing changes
    BOOST_CHECK(!onv.annihilateAll({1, 0}));
    BOOST_CHECK(onv_copy == onv);


    // Try to annihilate {5,2,0,1} and check if nothing changes
    BOOST_CHECK(!onv.annihilateAll({5, 2, 0, 1}));
    BOOST_CHECK(onv_copy == onv);
}


BOOST_AUTO_TEST_CASE ( annihilateAll_2 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)

    // Annihilate the indices {1, 2}
    GQCP::ONV ref_onv (6, 1, 16);  // "010000" (16)
    BOOST_CHECK(onv.annihilateAll({1, 2}));
    BOOST_CHECK(onv == ref_onv);
}


BOOST_AUTO_TEST_CASE ( annihilateAll_3 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)

    // Annihilate the indices {2, 1}
    GQCP::ONV ref_onv (6, 1, 16);  // "010000" (16)
    BOOST_CHECK(onv.annihilateAll({2, 1}));
    BOOST_CHECK(onv == ref_onv);
}


BOOST_AUTO_TEST_CASE ( annihilate_sign ) {

    // There should be a sign change when we annihilate on (lexical) index 2 for "10101" (21)
    GQCP::ONV onv1 (5,3, 21);
    int sign = 1;

    onv1.annihilate(2, sign);
    BOOST_CHECK_EQUAL(sign, -1);


    // There should be no sign change when we annihilate on (lexical) index 4 for "10101" (21)
    GQCP::ONV onv2 (5, 3, 21);
    sign = 1;

    onv2.annihilate(4, sign);
    BOOST_CHECK_EQUAL(sign, 1);

    // Annihilating the first occupied orbital on (lexical) index 0 for "10101" (21)
    GQCP::ONV onv3 (5, 3, 21);
    sign = 1;

    onv3.annihilate(0, sign);
    BOOST_CHECK_EQUAL(sign, 1);
}


BOOST_AUTO_TEST_CASE ( annihilateAll_sign_1 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)
    GQCP::ONV onv_copy = onv;  // to check if nothing happens
    int sign = 1;

    BOOST_CHECK_THROW(onv.annihilateAll({8, 4}, sign), std::invalid_argument);


    // Try to annihilate {0,1} and check if nothing changes
    BOOST_CHECK(!onv.annihilateAll({0, 1}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);


    // Try to annihilate {1,0} and check if nothing changes
    BOOST_CHECK(!onv.annihilateAll({1, 0}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);


    // Try to annihilate {5,2,0,1} and check if nothing changes
    BOOST_CHECK(!onv.annihilateAll({5, 2, 0, 1}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);
}


BOOST_AUTO_TEST_CASE ( annihilateAll_sign_2 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)

    // Annihilate the indices {1, 2}
    GQCP::ONV ref_onv (6, 1, 16);  // "010000" (16)
    int sign = 1;
    BOOST_CHECK(onv.annihilateAll({1, 2}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv == ref_onv);
}


BOOST_AUTO_TEST_CASE ( annihilateAll_sign_3 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)

    // Annihilate the indices {2, 1}
    GQCP::ONV ref_onv (6, 1, 16);  // "010000" (16)
    int sign = 1;
    BOOST_CHECK(onv.annihilateAll({2, 1}, sign));
    BOOST_CHECK(sign == -1);
    BOOST_CHECK(onv == ref_onv);
}


BOOST_AUTO_TEST_CASE ( create ) {

    GQCP::ONV onv (4, 1, 2);  // "0010" (2)

    // We shouldn't be able to create on index 9 (out of bounds)
    BOOST_CHECK_THROW(onv.create(9), std::invalid_argument);

    // We can't create on index 1
    BOOST_CHECK(!(onv.create(1)));

    // We can create on index 2
    BOOST_CHECK(onv.create(2));  // "0010" (2) -> "0110" (6)
    BOOST_CHECK_EQUAL(onv.get_unsigned_representation(), 6);
}


BOOST_AUTO_TEST_CASE ( createAll_1 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)
    GQCP::ONV onv_copy = onv;  // to check if nothing happens

    BOOST_CHECK_THROW(onv.createAll({8, 4}), std::invalid_argument);


    // Try to create {0,1} and check if nothing changes
    BOOST_CHECK(!onv.createAll({0, 1}));
    BOOST_CHECK(onv_copy == onv);


    // Try to create {5,2,0,1} and check if nothing changes
    BOOST_CHECK(!onv.createAll({5, 2, 0, 1}));
    BOOST_CHECK(onv_copy == onv);
}


BOOST_AUTO_TEST_CASE ( createAll_2 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)

    // Create the indices {0, 3}
    GQCP::ONV ref_onv (6, 5, 31);  // "011111" (31)
    BOOST_CHECK(onv.createAll({0, 3}));
    BOOST_CHECK(onv == ref_onv);
}


BOOST_AUTO_TEST_CASE ( create_sign ) {

    // There should be a sign change when we create on (lexical) index 1 for "10101" (21)
    GQCP::ONV onv1 (5, 3, 21);
    int sign = 1;

    onv1.create(1, sign);
    BOOST_CHECK_EQUAL(sign, -1);


    // There should be no sign change when we create on (lexical) index 3 for "10101" (21)
    GQCP::ONV onv2 (5, 3, 21);
    sign = 1;

    onv2.create(3, sign);
    BOOST_CHECK_EQUAL(sign, 1);

    // Creating the last orbital on (lexical) index 4 for "00101" (5)
    GQCP::ONV onv3 (5, 2, 5);
    sign = 1;

    onv3.create(4, sign);
    BOOST_CHECK_EQUAL(sign, 1);

    // Creating the first orbital on (lexical) index 0 for "10100" (5)
    GQCP::ONV onv4 (5, 2, 20);
    sign = 1;

    onv4.create(0, sign);
    BOOST_CHECK_EQUAL(sign, 1);
}


BOOST_AUTO_TEST_CASE ( createAll_sign_1 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)
    GQCP::ONV onv_copy = onv;  // to check if nothing happens

    BOOST_CHECK_THROW(onv.createAll({8, 4}), std::invalid_argument);


    // Try to create {0,1} and check if nothing changes
    int sign = 1;
    BOOST_CHECK(!onv.createAll({0, 1}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);


    // Try to create {5,2,0,1} and check if nothing changes
    BOOST_CHECK(!onv.createAll({5, 2, 0, 1}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv_copy == onv);
}


BOOST_AUTO_TEST_CASE ( createAll_sign_2 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)

    // Create the indices {0, 3}
    GQCP::ONV ref_onv (6, 5, 31);  // "011111" (31)
    int sign = 1;
    BOOST_CHECK(onv.createAll({0, 3}, sign));
    BOOST_CHECK(sign == -1);
    BOOST_CHECK(onv == ref_onv);
}


BOOST_AUTO_TEST_CASE ( createAll_sign_3 ) {

    GQCP::ONV onv (6, 3, 22);  // "010110" (22)

    // Create the indices {3, 0}
    int sign = 1;
    GQCP::ONV ref_onv (6, 5, 31);  // "011111" (31)
    BOOST_CHECK(onv.createAll({3, 0}, sign));
    BOOST_CHECK(sign == 1);
    BOOST_CHECK(onv == ref_onv);
}


BOOST_AUTO_TEST_CASE ( countNumberOfDifferences ) {

    GQCP::ONV onv1 (5, 3, 21);  // "10101" (21)
    GQCP::ONV onv2 (5, 3, 22);  // "10110" (22)
    GQCP::ONV onv3 (5, 3, 26);  // "11010" (26)

    BOOST_CHECK_EQUAL(onv1.countNumberOfDifferences(onv1), 0);
    BOOST_CHECK_EQUAL(onv2.countNumberOfDifferences(onv2), 0);
    BOOST_CHECK_EQUAL(onv3.countNumberOfDifferences(onv3), 0);

    BOOST_CHECK_EQUAL(onv1.countNumberOfDifferences(onv2), 2);
    BOOST_CHECK_EQUAL(onv1.countNumberOfDifferences(onv3), 4);
    BOOST_CHECK_EQUAL(onv2.countNumberOfDifferences(onv3), 2);
}


BOOST_AUTO_TEST_CASE ( findDifferentOccupations ) {

    GQCP::ONV onv1 (5, 3, 21);  // "10101" (21)
    GQCP::ONV onv2 (5, 3, 22);  // "10110" (22)
    GQCP::ONV onv3 (5, 3, 26);  // "11010" (26)

    BOOST_TEST(onv1.findDifferentOccupations(onv2) == (std::vector<size_t> {0}), boost::test_tools::per_element());
    BOOST_TEST(onv2.findDifferentOccupations(onv1) == (std::vector<size_t> {1}), boost::test_tools::per_element());

    BOOST_TEST(onv1.findDifferentOccupations(onv3) == (std::vector<size_t> {0, 2}), boost::test_tools::per_element());
    BOOST_TEST(onv3.findDifferentOccupations(onv1) == (std::vector<size_t> {1, 3}), boost::test_tools::per_element());

    BOOST_TEST(onv2.findDifferentOccupations(onv3) == (std::vector<size_t> {2}), boost::test_tools::per_element());
    BOOST_TEST(onv3.findDifferentOccupations(onv2) == (std::vector<size_t> {3}), boost::test_tools::per_element());
}


BOOST_AUTO_TEST_CASE ( findMatchingOccupations ) {

    GQCP::ONV onv1 (5, 3, 21);  // "10101" (21)
    GQCP::ONV onv2 (5, 3, 22);  // "10110" (22)
    GQCP::ONV onv3 (5, 3, 26);  // "11010" (26)

    BOOST_TEST(onv1.findMatchingOccupations(onv2) == (std::vector<size_t> {2,4}), boost::test_tools::per_element());
    BOOST_TEST(onv2.findMatchingOccupations(onv1) == (std::vector<size_t> {2,4}), boost::test_tools::per_element());

    BOOST_TEST(onv1.findMatchingOccupations(onv3) == (std::vector<size_t> {4}), boost::test_tools::per_element());
    BOOST_TEST(onv3.findMatchingOccupations(onv1) == (std::vector<size_t> {4}), boost::test_tools::per_element());

    BOOST_TEST(onv2.findMatchingOccupations(onv3) == (std::vector<size_t> {1,4}), boost::test_tools::per_element());
    BOOST_TEST(onv3.findMatchingOccupations(onv2) == (std::vector<size_t> {1,4}), boost::test_tools::per_element());
}
