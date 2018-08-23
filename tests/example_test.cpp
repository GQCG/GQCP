#define BOOST_TEST_MODULE "example_test"


#include "temp.hpp"



#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( isTrue_true ) {

    BOOST_CHECK(true);
}
