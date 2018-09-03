#define BOOST_TEST_MODULE "Atom"


#include "Atom.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain



BOOST_AUTO_TEST_CASE ( Atom_constructor ) {

    GQCG::Atom atom {1, 0.0, 0.1, 0.2};
}


BOOST_AUTO_TEST_CASE ( Atom_isSmallerThan ) {

    GQCG::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom2 {2, 0.0, 0.1, 0.2};
    GQCG::Atom atom3 {2, 0.1, 0.2, 0.2};
    GQCG::Atom atom4 {2, 0.1, 0.2, 0.3};


    // Check if operator< does what is expected
    BOOST_CHECK(atom1.isSmallerThan(atom2));
    BOOST_CHECK(atom2.isSmallerThan(atom3));
    BOOST_CHECK(atom3.isSmallerThan(atom4));

    BOOST_CHECK(!(atom2.isSmallerThan(atom1)));
    BOOST_CHECK(!(atom3.isSmallerThan(atom2)));
    BOOST_CHECK(!(atom4.isSmallerThan(atom3)));


    // Check if the tolerance works
    BOOST_CHECK(!(atom2.isSmallerThan(atom3, 0.2)));
    BOOST_CHECK(!(atom3.isSmallerThan(atom2, 0.2)));
}


BOOST_AUTO_TEST_CASE ( Atom_operator_smaller_than ) {

    GQCG::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom2 {2, 0.0, 0.1, 0.2};

    // A small test to check if we can operator<
    BOOST_CHECK(atom1 < atom2);
}


BOOST_AUTO_TEST_CASE ( Atom_isEqualTo ) {

    GQCG::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom2 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom3 {2, 0.0, 0.1, 0.2};
    GQCG::Atom atom4 {1, 0.1, 0.2, 0.3};

    // Check if they're equal
    BOOST_CHECK(atom1.isEqualTo(atom2));

    // Check if different atomic numbers cause inequality
    BOOST_CHECK(!(atom1.isEqualTo(atom3)));

    // Check if different coordinates cause inequality
    BOOST_CHECK(!(atom1.isEqualTo(atom4)));


    // Check if the tolerance works
    BOOST_CHECK(atom1.isEqualTo(atom2, 0.2));
}


BOOST_AUTO_TEST_CASE ( Atom_operator_equals ) {

    GQCG::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom2 {1, 0.0, 0.1, 0.2};

    // A small test to check if we can operator==
    BOOST_CHECK(atom1 == atom2);
}


BOOST_AUTO_TEST_CASE ( Atom_operator_ostream ) {

    GQCG::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom2 {2, 0.1, 0.2, 0.3};


    std::cout << atom1 << std::endl;
    std::cout << atom2 << std::endl;
}


BOOST_AUTO_TEST_CASE ( calculateDistance ) {

    // Create some atoms
    GQCG::Atom atom1 {1, 0, 3, 0};
    GQCG::Atom atom2 {1, 0, 0, 4};
    GQCG::Atom atom3 {1, 3, 0, 0};
    GQCG::Atom atom4 {1, 0, 0, 5};


    // Check their distances
    BOOST_CHECK(std::abs(atom1.calculateDistance(atom2) - 5) < 1.0e-12);
    BOOST_CHECK(std::abs(atom1.calculateDistance(atom3) - std::sqrt(18.0)) < 1.0e-12);
    BOOST_CHECK(std::abs(atom1.calculateDistance(atom2) - atom2.calculateDistance(atom3)) < 1.0e-12);
    BOOST_CHECK(std::abs(atom2.calculateDistance(atom3) - 5) < 1.0e-12);
    BOOST_CHECK(std::abs(atom2.calculateDistance(atom4) - 1) < 1.0e-12);

    // Check that the distances are symmetric
    BOOST_CHECK(std::abs(atom1.calculateDistance(atom2) - atom2.calculateDistance(atom1)) < 1.0e-12);
    BOOST_CHECK(std::abs(atom1.calculateDistance(atom3) - atom3.calculateDistance(atom1)) < 1.0e-12);
    BOOST_CHECK(std::abs(atom2.calculateDistance(atom3) - atom3.calculateDistance(atom2)) < 1.0e-12);
}
