#define BOOST_TEST_MODULE "FockSpaceProduct"


#include "FockSpace/FockSpaceProduct.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( FockSpaceProduct_constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::FockSpaceProduct (10, 5, 5));
}


BOOST_AUTO_TEST_CASE ( FockSpaceProduct_dimension) {

    BOOST_CHECK_EQUAL(GQCP::FockSpaceProduct::calculateDimension(10, 1, 1), 100);
    BOOST_CHECK_EQUAL(GQCP::FockSpaceProduct::calculateDimension(6, 2, 2), 225);
    BOOST_CHECK_EQUAL(GQCP::FockSpaceProduct::calculateDimension(8, 3, 3), 3136);

    BOOST_CHECK_EQUAL(GQCP::FockSpaceProduct::calculateDimension(10, 2, 0), 45);
    BOOST_CHECK_EQUAL(GQCP::FockSpaceProduct::calculateDimension(6, 3, 1), 120);
    BOOST_CHECK_EQUAL(GQCP::FockSpaceProduct::calculateDimension(8, 4, 2), 1960);
}
