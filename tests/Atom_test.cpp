#define BOOST_TEST_MODULE "Atom"


#include "Atom.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain



BOOST_AUTO_TEST_CASE ( Atom_constructor ) {

    GQCG::Atom atom {1, 0.0, 0.1, 0.2};
}


BOOST_AUTO_TEST_CASE ( Atom_operator_equals ) {

    GQCG::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom2 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom3 {2, 0.0, 0.1, 0.2};
    GQCG::Atom atom4 {1, 0.1, 0.2, 0.3};

    // Check if they're equal
    BOOST_CHECK(atom1 == atom2);

    // Check if different atomic numbers cause inequality
    BOOST_CHECK(!(atom1 == atom3));

    // Check if different coordinates cause inequality
    BOOST_CHECK(!(atom1 == atom4));
}


BOOST_AUTO_TEST_CASE ( Atom_operator_smaller_than ) {

    GQCG::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom2 {2, 0.0, 0.1, 0.2};
    GQCG::Atom atom3 {2, 0.1, 0.2, 0.2};
    GQCG::Atom atom4 {2, 0.1, 0.2, 0.3};


    // Check if operator< does what is expected
    BOOST_CHECK(atom1 < atom2);
    BOOST_CHECK(atom2 < atom3);
    BOOST_CHECK(atom3 < atom4);

    BOOST_CHECK(!(atom2 < atom1));
    BOOST_CHECK(!(atom3 < atom2));
    BOOST_CHECK(!(atom4 < atom3));
}


BOOST_AUTO_TEST_CASE ( Atom_operator_ostream ) {

    GQCG::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom2 {2, 0.1, 0.2, 0.3};


    std::cout << atom1 << std::endl;
    std::cout << atom2 << std::endl;
}
