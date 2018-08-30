#define BOOST_TEST_MODULE "Atom"


#include "Atom.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain



BOOST_AUTO_TEST_CASE ( constructor ) {

    GQCG::Atom atom {1, 0.0, 0.1, 0.2};
}


BOOST_AUTO_TEST_CASE ( Atom_isEqualTo ) {

    GQCG::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom2 {1, 0.0, 0.1, 0.2};
    GQCG::Atom atom3 {2, 0.0, 0.1, 0.2};
    GQCG::Atom atom4 {1, 0.1, 0.2, 0.3};

    // Check if they're equal
    BOOST_CHECK(atom1.isEqualTo(atom2));

    // Check if different atomic numbers cause inequality
    BOOST_CHECK(!atom1.isEqualTo(atom3));

    // Check if different coordinates cause inequality
    BOOST_CHECK(!atom1.isEqualTo(atom4));

    // Check if the tolerance works as expected
    BOOST_CHECK(atom1.isEqualTo(atom4, 0.5));
}
