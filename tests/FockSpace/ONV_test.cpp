// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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


#include "FockSpace/ONV.hpp"
#include "FockSpace/FockSpace.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE ( operator_equals_onv ) {

    GQCP::ONV spin_string1 (6, 3, 19);  // "010011" (19)
    GQCP::ONV spin_string2 (6, 3, 19);  // "010011" (19)
    GQCP::ONV spin_string3 (7, 3, 19);  // "0010011" (19)

    BOOST_CHECK(spin_string1 == spin_string2);
    BOOST_CHECK(!(spin_string1 == spin_string3));  // wrong number of orbitals
}


BOOST_AUTO_TEST_CASE ( ONV_constructor ) {
    // Create ONV : 10 considered bits and 5 set bits with distributed as "0000011111" = 31
    GQCP::ONV onv1 (10, 5, 31);

    // Create Fock space with 10 orbitals and 5 electrons
    GQCP::FockSpace fock_space (10, 5);

    // Ask for the first ONV in reverse lexicographic order : "0000011111" = 31
    GQCP::ONV onv2 = fock_space.get_ONV(0);

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


BOOST_AUTO_TEST_CASE (  operator_not_equals_onv ) {

    GQCP::ONV spin_string1 (6, 3, 19);  // "010011" (19)
    GQCP::ONV spin_string2 (6, 3, 19);  // "010011" (19)
    GQCP::ONV spin_string3 (7, 3, 19);  // "0010011" (19)

    BOOST_CHECK(!(spin_string1 != spin_string2));
    BOOST_CHECK(spin_string1 != spin_string3);  // wrong number of orbitals
}


BOOST_AUTO_TEST_CASE ( address_onv ) {

    GQCP::FockSpace fock_space (6, 3);

    // The address of the string "010011" (19) should be 4
    GQCP::ONV onv (6, 3, 19);

    BOOST_CHECK_EQUAL(fock_space.getAddress(onv), 4);
}


BOOST_AUTO_TEST_CASE ( set_Next_onv ) {
    GQCP::FockSpace fock_space(5, 3);
    // K = 5, N = 3 <-> "00111"
    GQCP::ONV spin_string = fock_space.get_ONV(0);
    // The lexical permutations are: "00111" (7), "01011" (11), "01101" (13), "01110" (14), etc.

    // Check permutations one after the other

    fock_space.setNext(spin_string);  // "01011" (11)
    BOOST_CHECK_EQUAL(spin_string.get_unsigned_representation(), 11);
    GQCP::VectorXs x1 (3);
    x1 << 0, 1, 3;
    BOOST_CHECK(x1.isApprox(spin_string.get_occupation_indices()));

    fock_space.setNext(spin_string);  // "01101" (13)
    BOOST_CHECK_EQUAL(spin_string.get_unsigned_representation(), 13);
    GQCP::VectorXs x2 (3);
    x2 << 0, 2, 3;
    BOOST_CHECK(x2.isApprox(spin_string.get_occupation_indices()));

    fock_space.setNext(spin_string);  // "01110" (14)
    BOOST_CHECK_EQUAL(spin_string.get_unsigned_representation(), 14);
    GQCP::VectorXs x3 (3);
    x3 << 1, 2, 3;
    BOOST_CHECK(x3.isApprox(spin_string.get_occupation_indices()));
}


BOOST_AUTO_TEST_CASE ( slice_onv ) {

    // Check for throws
    GQCP::ONV spin_string (4, 2, 5);  // "0101" (5)
    BOOST_CHECK_THROW(spin_string.slice(2, 1), std::invalid_argument);  // index_end should be larger than index_start
    BOOST_CHECK_THROW(spin_string.slice(2, 2), std::invalid_argument);  // index_end should be larger than index_start
    BOOST_CHECK_THROW(spin_string.slice(2, 6), std::invalid_argument);  // index_end is out of bounds


    GQCP::ONV spin_string1 (7, 3, 41);
    BOOST_CHECK_EQUAL(spin_string1.slice(2, 6), 10);  // "0[1010]01" (41) -> "1010" (10)

    GQCP::ONV spin_string2 (7, 5, 115);
    BOOST_CHECK_EQUAL(spin_string2.slice(2, 7), 28);  // "[11100]11" (115) -> "11100" (28)

    GQCP::ONV spin_string3 (7, 5, 115);
    BOOST_CHECK_EQUAL(spin_string3.slice(2, 5), 4);  // "11[100]11" (115) -> "100" (4)

    GQCP::ONV spin_string4 (6, 3, 19);
    BOOST_CHECK_EQUAL(spin_string4.slice(2, 6), 4);  // "[0100]11" (18) -> "0100" (4)

    GQCP::ONV spin_string5 (6, 3, 19);
    BOOST_CHECK_EQUAL(spin_string5.slice(5, 6), 0);  // "[0]10011" (18) -> "0" (0)

    GQCP::ONV spin_string6 (6, 3, 19);
    BOOST_CHECK_EQUAL(spin_string6.slice(4, 5), 1);  // "0[1]0011" (19) -> "1" (1)
}


BOOST_AUTO_TEST_CASE ( isOccupied_onv ) {

    GQCP::ONV spin_string (4, 1, 2);  // "0010" (2)

    // We shouldn't be able to check on index 9 (out of bounds)
    BOOST_CHECK_THROW(spin_string.isOccupied(4), std::invalid_argument);

    // Index 1 is occupied in "0010" (2)
    BOOST_CHECK(spin_string.isOccupied(1));

    // Index 2 is not occupied "0010" (2)
    BOOST_CHECK(!spin_string.isOccupied(2));
}


BOOST_AUTO_TEST_CASE ( annihilate_onv ) {

    GQCP::ONV spin_string (4, 2, 10);  // "1010" (10)

    // We shouldn't be able to annihilate on index 5 (out of bounds)
    BOOST_CHECK_THROW(spin_string.annihilate(5), std::invalid_argument);

    // We can't annihilate on index 2
    BOOST_CHECK(!(spin_string.annihilate(2)));

    // We can annihilate on index 1
    BOOST_CHECK(spin_string.annihilate(1));  // "1010" (10) -> "1000" (8)
    BOOST_CHECK_EQUAL(spin_string.get_unsigned_representation(), 8);
    // Test if updating throws an error (no longer 2 electrons)
    BOOST_CHECK_THROW(spin_string.updateOccupationIndices(), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( annihilate_sign_onv ) {

    // There should be a sign change when we annihilate on (lexical) index 2 for "10101" (21)
    GQCP::ONV spin_string1 (5,3, 21);
    int sign = 1;

    spin_string1.annihilate(2, sign);
    BOOST_CHECK_EQUAL(sign, -1);


    // There should be no sign change when we annihilate on (lexical) index 4 for "10101" (21)
    GQCP::ONV spin_string2 (5, 3, 21);
    sign = 1;

    spin_string2.annihilate(4, sign);
    BOOST_CHECK_EQUAL(sign, 1);

    // Annihilating the first occupied SO on (lexical) index 0 for "10101" (21)
    GQCP::ONV spin_string3 (5, 3, 21);
    sign = 1;

    spin_string3.annihilate(0, sign);
    BOOST_CHECK_EQUAL(sign, 1);
}


BOOST_AUTO_TEST_CASE ( create_onv ) {

    GQCP::ONV spin_string (4, 1, 2);  // "0010" (2)

    // We shouldn't be able to create on index 9 (out of bounds)
    BOOST_CHECK_THROW(spin_string.create(9), std::invalid_argument);

    // We can't create on index 1
    BOOST_CHECK(!(spin_string.create(1)));

    // We can create on index 2
    BOOST_CHECK(spin_string.create(2));  // "0010" (2) -> "0110" (6)
    BOOST_CHECK_EQUAL(spin_string.get_unsigned_representation(), 6);
}


BOOST_AUTO_TEST_CASE ( create_sign_onv ) {

    // There should be a sign change when we create on (lexical) index 1 for "10101" (21)
    GQCP::ONV spin_string1 (5, 3, 21);
    int sign = 1;

    spin_string1.create(1, sign);
    BOOST_CHECK_EQUAL(sign, -1);


    // There should be no sign change when we create on (lexical) index 3 for "10101" (21)
    GQCP::ONV spin_string2 (5, 3, 21);
    sign = 1;

    spin_string2.create(3, sign);
    BOOST_CHECK_EQUAL(sign, 1);

    // Creating the last SO on (lexical) index 4 for "00101" (5)
    GQCP::ONV spin_string3 (5, 2, 5);
    sign = 1;

    spin_string3.create(4, sign);
    BOOST_CHECK_EQUAL(sign, 1);

    // Creating the first SO on (lexical) index 0 for "10100" (5)
    GQCP::ONV spin_string4 (5, 2, 20);
    sign = 1;

    spin_string4.create(0, sign);
    BOOST_CHECK_EQUAL(sign, 1);
}


BOOST_AUTO_TEST_CASE ( ONV_address_setNext_fullspace ) {
    // Here we will test a full permutation through a Fock space of K = 15, N = 5
    GQCP::FockSpace fock_space (15, 5);

    // Retrieve the first ONV of the Fock space
    GQCP::ONV onv_test = fock_space.get_ONV(0);

    const size_t dimension_fock_space = 3003;
    bool is_correct = true;  // variable that is updated to false if an unexpected result occurs

    // Iterate through the Fock space in reverse lexicographical order and test whether address matches
    for (int i = 0; i<dimension_fock_space; i++) {

        // Tests address
        if (i != fock_space.getAddress(onv_test)) {
            is_correct = false;
        }

        // transforms the given ONV to the next ONV in the Fock space
        if (i < dimension_fock_space - 1) {
            fock_space.setNext(onv_test);
        }
    }

    // Checks if no unexpected results occured in a full iteration
    BOOST_CHECK(is_correct);
}


BOOST_AUTO_TEST_CASE ( countNumberOfDifferences ) {

    GQCP::ONV spin_string1 (5, 3, 21);  // "10101" (21)
    GQCP::ONV spin_string2 (5, 3, 22);  // "10110" (22)
    GQCP::ONV spin_string3 (5, 3, 26);  // "11010" (26)

    BOOST_CHECK_EQUAL(spin_string1.countNumberOfDifferences(spin_string1), 0);
    BOOST_CHECK_EQUAL(spin_string2.countNumberOfDifferences(spin_string2), 0);
    BOOST_CHECK_EQUAL(spin_string3.countNumberOfDifferences(spin_string3), 0);

    BOOST_CHECK_EQUAL(spin_string1.countNumberOfDifferences(spin_string2), 2);
    BOOST_CHECK_EQUAL(spin_string1.countNumberOfDifferences(spin_string3), 4);
    BOOST_CHECK_EQUAL(spin_string2.countNumberOfDifferences(spin_string3), 2);
}


BOOST_AUTO_TEST_CASE ( findOccupiedDifferences ) {

    GQCP::ONV spin_string1 (5, 3, 21);  // "10101" (21)
    GQCP::ONV spin_string2 (5, 3, 22);  // "10110" (22)
    GQCP::ONV spin_string3 (5, 3, 26);  // "11010" (26)

    BOOST_TEST(spin_string1.findDifferentOccupations(spin_string2) == (std::vector<size_t> {0}), boost::test_tools::per_element());
    BOOST_TEST(spin_string2.findDifferentOccupations(spin_string1) == (std::vector<size_t> {1}), boost::test_tools::per_element());

    BOOST_TEST(spin_string1.findDifferentOccupations(spin_string3) == (std::vector<size_t> {0, 2}), boost::test_tools::per_element());
    BOOST_TEST(spin_string3.findDifferentOccupations(spin_string1) == (std::vector<size_t> {1, 3}), boost::test_tools::per_element());

    BOOST_TEST(spin_string2.findDifferentOccupations(spin_string3) == (std::vector<size_t> {2}), boost::test_tools::per_element());
    BOOST_TEST(spin_string3.findDifferentOccupations(spin_string2) == (std::vector<size_t> {3}), boost::test_tools::per_element());
}


BOOST_AUTO_TEST_CASE ( findMatchingOccupations ) {

    GQCP::ONV spin_string1 (5, 3, 21);  // "10101" (21)
    GQCP::ONV spin_string2 (5, 3, 22);  // "10110" (22)
    GQCP::ONV spin_string3 (5, 3, 26);  // "11010" (26)

    BOOST_TEST(spin_string1.findMatchingOccupations(spin_string2) == (std::vector<size_t> {2,4}), boost::test_tools::per_element());
    BOOST_TEST(spin_string2.findMatchingOccupations(spin_string1) == (std::vector<size_t> {2,4}), boost::test_tools::per_element());

    BOOST_TEST(spin_string1.findMatchingOccupations(spin_string3) == (std::vector<size_t> {4}), boost::test_tools::per_element());
    BOOST_TEST(spin_string3.findMatchingOccupations(spin_string1) == (std::vector<size_t> {4}), boost::test_tools::per_element());

    BOOST_TEST(spin_string2.findMatchingOccupations(spin_string3) == (std::vector<size_t> {1,4}), boost::test_tools::per_element());
    BOOST_TEST(spin_string3.findMatchingOccupations(spin_string2) == (std::vector<size_t> {1,4}), boost::test_tools::per_element());
}