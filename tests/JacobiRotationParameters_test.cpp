#define BOOST_TEST_MODULE "JacobiRotationParameters"


#include "JacobiRotationParameters.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( JacobiRotationParameters_constructor ) {

    // Check if a correct constructor works
    GQCP::JacobiRotationParameters jacobi_parameters (3, 1, 0.5);

    // Check if we can't construct when p < q
    BOOST_CHECK_THROW(GQCP::JacobiRotationParameters(1, 3, 0.5), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( operator_stream ) {

    GQCP::JacobiRotationParameters jacobi_parameters (3, 1, 0.5);
    std::cout << jacobi_parameters << std::endl;
}
