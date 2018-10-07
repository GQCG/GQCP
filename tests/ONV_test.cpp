#define BOOST_TEST_MODULE "ONV"


#include "ONV.hpp"
#include "FockSpace/FockSpace.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain



BOOST_AUTO_TEST_CASE ( ONV_constructor ) {
    // Create ONV with 10 considered bits and 5 set bits with distributed as "0000011111" = 31
    GQCG::ONV onv1 (10,5,31);

    // Create Fock Space with 10 orbitals and 5 electrons
    GQCG::FockSpace fock_space(10,5);

    // Ask the first ONV in reverse lexicographic order : "0000011111" = 31
    GQCG::ONV onv2 = fock_space.get_ONV(0);

    // Check if both unsigned representations match (are equal to 31 and have the same considered bits)
    BOOST_CHECK(onv1 == onv2);
    BOOST_CHECK(onv1.get_urepresentation() == 31);

    // Test if the occupation numbers are correct.
    GQCG::VectorXs x(5);
    x << 0, 1, 2, 3, 4;
    BOOST_CHECK(x.isApprox(onv1.get_occupations()));

    // Test if setting incompatible representation throws an error
    BOOST_CHECK_THROW(onv1.set_representation(1), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( operator_equals_onv ) {

    GQCG::ONV spin_string1 (6, 3, 19);  // "010011" (19)
    GQCG::ONV spin_string2 (6, 3, 19);  // "010011" (19)
    GQCG::ONV spin_string3 (7, 3, 19);  // "0010011" (19)

    BOOST_CHECK(spin_string1 == spin_string2);
    BOOST_CHECK(!(spin_string1 == spin_string3));  // wrong number of orbitals
}

BOOST_AUTO_TEST_CASE (  operator_not_equals_onv ) {

    GQCG::ONV spin_string1 (6, 3, 19);  // "010011" (19)
    GQCG::ONV spin_string2 (6, 3, 19);  // "010011" (19)
    GQCG::ONV spin_string3 (7, 3, 19);  // "0010011" (19)

    BOOST_CHECK(!(spin_string1 != spin_string2));
    BOOST_CHECK(spin_string1 != spin_string3);  // wrong number of orbitals
}

BOOST_AUTO_TEST_CASE ( address_onv ) {

    GQCG::FockSpace fock_space(6,3);

    // The address of the string "010011" (19) should be 4
    GQCG::ONV onv(6,3,19);

    BOOST_CHECK_EQUAL(fock_space.getAddress(onv), 4);
}

BOOST_AUTO_TEST_CASE ( set_Next_onv ) {
    GQCG::FockSpace fock_space(5,3);
    // K = 5, N = 3 <-> "00111"
    GQCG::ONV spin_string = fock_space.get_ONV(0);
    // The lexical permutations are: "00111" (7), "01011" (11), "01101" (13), "01110" (14), etc.

    // Check permutations one after the other

    fock_space.setNext(spin_string);  // "01011" (11)
    BOOST_CHECK_EQUAL(spin_string.get_urepresentation(), 11);
    GQCG::VectorXs x1(3);
    x1 << 0, 1, 3;
    BOOST_CHECK(x1.isApprox(spin_string.get_occupations()));

    fock_space.setNext(spin_string);  // "01101" (13)
    BOOST_CHECK_EQUAL(spin_string.get_urepresentation(), 13);
    GQCG::VectorXs x2(3);
    x2 << 0, 2, 3;
    BOOST_CHECK(x2.isApprox(spin_string.get_occupations()));

    fock_space.setNext(spin_string);  // "01110" (14)
    BOOST_CHECK_EQUAL(spin_string.get_urepresentation(), 14);
    GQCG::VectorXs x3(3);
    x3 << 1, 2, 3;
    BOOST_CHECK(x3.isApprox(spin_string.get_occupations()));
}


BOOST_AUTO_TEST_CASE ( slice_onv ) {

    // Check for throws
    GQCG::ONV spin_string (4,2,5);  // "0101" (5)
    BOOST_CHECK_THROW(spin_string.slice(2, 1), std::invalid_argument);  // index_end should be larger than index_start
    BOOST_CHECK_THROW(spin_string.slice(2, 2), std::invalid_argument);  // index_end should be larger than index_start
    BOOST_CHECK_THROW(spin_string.slice(2, 6), std::invalid_argument);  // index_end is out of bounds


    GQCG::ONV spin_string1 (7,3,41);
    BOOST_CHECK_EQUAL(spin_string1.slice(2, 6), 10);  // "0[1010]01" (41) -> "1010" (10)

    GQCG::ONV spin_string2 (7,5,115);
    BOOST_CHECK_EQUAL(spin_string2.slice(2, 7), 28);  // "[11100]11" (115) -> "11100" (28)

    GQCG::ONV spin_string3 (7,5, 115);
    BOOST_CHECK_EQUAL(spin_string3.slice(2, 5), 4);  // "11[100]11" (115) -> "100" (4)

    GQCG::ONV spin_string4 (6,3,19);
    BOOST_CHECK_EQUAL(spin_string4.slice(2, 6), 4);  // "[0100]11" (18) -> "0100" (4)

    GQCG::ONV spin_string5 (6,3,19);
    BOOST_CHECK_EQUAL(spin_string5.slice(5, 6), 0);  // "[0]10011" (18) -> "0" (0)

    GQCG::ONV spin_string6 (6,3, 19);
    BOOST_CHECK_EQUAL(spin_string6.slice(4, 5), 1);  // "0[1]0011" (19) -> "1" (1)
}


BOOST_AUTO_TEST_CASE ( isOccupied_onv ) {

    GQCG::ONV spin_string (4,1, 2);  // "0010" (2)

    // We shouldn't be able to check on index 9 (out of bounds)
    BOOST_CHECK_THROW(spin_string.isOccupied(4), std::invalid_argument);

    // Index 1 is occupied in "0010" (2)
    BOOST_CHECK(spin_string.isOccupied(1));

    // Index 2 is not occupied "0010" (2)
    BOOST_CHECK(!spin_string.isOccupied(2));
}


BOOST_AUTO_TEST_CASE ( annihilate_onv ) {

    GQCG::ONV spin_string (4,2,10);  // "1010" (10)

    // We shouldn't be able to annihilate on index 5 (out of bounds)
    BOOST_CHECK_THROW(spin_string.annihilate(5), std::invalid_argument);

    // We can't annihilate on index 2
    BOOST_CHECK(!(spin_string.annihilate(2)));

    // We can annihilate on index 1
    BOOST_CHECK(spin_string.annihilate(1));  // "1010" (10) -> "1000" (8)
    BOOST_CHECK_EQUAL(spin_string.get_urepresentation(), 8);
    // Test if updating throws an error (no longer 2 electrons)
    BOOST_CHECK_THROW(spin_string.update(), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( annihilate_sign_onv ) {

    // There should be a sign change when we annihilate on (lexical) index 2 for "10101" (21)
    GQCG::ONV spin_string1 (5,3, 21);
    int sign = 1;

    spin_string1.annihilate(2, sign);
    BOOST_CHECK_EQUAL(sign, -1);


    // There should be no sign change when we annihilate on (lexical) index 4 for "10101" (21)
    GQCG::ONV spin_string2 (5,3, 21);
    sign = 1;

    spin_string2.annihilate(4, sign);
    BOOST_CHECK_EQUAL(sign, 1);

    // Boundary condition, annihilating the first occupied SO on (lexical) index 0 for "10101" (21)
    GQCG::ONV spin_string3 (5,3, 21);
    sign = 1;

    spin_string3.annihilate(0, sign);
    BOOST_CHECK_EQUAL(sign, 1);
}


BOOST_AUTO_TEST_CASE ( create_onv ) {

    GQCG::ONV spin_string (4,1,2);  // "0010" (2)

    // We shouldn't be able to create on index 9 (out of bounds)
    BOOST_CHECK_THROW(spin_string.create(9), std::invalid_argument);

    // We can't create on index 1
    BOOST_CHECK(!(spin_string.create(1)));

    // We can create on index 2
    BOOST_CHECK(spin_string.create(2));  // "0010" (2) -> "0110" (6)
    BOOST_CHECK_EQUAL(spin_string.get_urepresentation(), 6);
}


BOOST_AUTO_TEST_CASE ( create_sign_onv ) {

    // There should be a sign change when we create on (lexical) index 1 for "10101" (21)
    GQCG::ONV spin_string1 (5,3, 21);
    int sign = 1;

    spin_string1.create(1, sign);
    BOOST_CHECK_EQUAL(sign, -1);


    // There should be no sign change when we create on (lexical) index 3 for "10101" (21)
    GQCG::ONV spin_string2 (5,3, 21);
    sign = 1;

    spin_string2.create(3, sign);
    BOOST_CHECK_EQUAL(sign, 1);

    // Boundary condition, creating the last SO on (lexical) index 4 for "00101" (5)
    GQCG::ONV spin_string3 (5,2, 5);
    sign = 1;

    spin_string3.create(4, sign);
    BOOST_CHECK_EQUAL(sign, 1);

    // Boundary condition, creating the first SO on (lexical) index 0 for "10100" (5)
    GQCG::ONV spin_string4 (5,2,20);
    sign = 1;

    spin_string4.create(0, sign);
    BOOST_CHECK_EQUAL(sign, 1);
}


BOOST_AUTO_TEST_CASE ( ONV_address_setNext_fullspace ) {
    // Here we will test a full permutation through a Fock space of K = 15, N = 5
    GQCG::FockSpace fock_space(15,5);

    // Retrieve the first ONV of the Fock space
    GQCG::ONV onv_test = fock_space.get_ONV(0);

    const size_t dimension_fock_space = 3003;
    bool is_correct = true;  // variable that is updated to false if an unexpected result occurs

    // Iterate through the Fock space in reverse lexicographical order and test whether address matches
    for(int i = 0; i<dimension_fock_space; i++) {

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
