#define BOOST_TEST_MODULE "AOBasis"


#include "AOBasis.hpp"

#include "Molecule.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( AOBasis_constructor ) {

    // Check if we can construct an AOBasis object
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    GQCG::AOBasis basis (water, "STO-3G");
}


BOOST_AUTO_TEST_CASE ( number_of_basis_functions ) {

    // Check the number of basis functions in water
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    GQCG::AOBasis basis (water, "STO-3G");

    BOOST_CHECK_EQUAL(basis.get_number_of_basis_functions(), 7);
}
