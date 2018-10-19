#define BOOST_TEST_MODULE "elements"


#include "elements.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( elementToAtomicNumber ) {
    
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("H"), 1);
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("C"), 6);
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("Mg"), 12);
}


BOOST_AUTO_TEST_CASE ( atomicNumberToElement ) {

    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(1), "H");
    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(6), "C");
    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(12), "Mg");
}
