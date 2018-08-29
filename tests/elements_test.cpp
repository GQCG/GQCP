#define BOOST_TEST_MODULE "elements"


#include "elements.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( element_to_atomic_number ) {
    
    BOOST_CHECK_EQUAL(GQCG::elements::element_to_atomic_number("H"), 1);
    BOOST_CHECK_EQUAL(GQCG::elements::element_to_atomic_number("C"), 6);

}
