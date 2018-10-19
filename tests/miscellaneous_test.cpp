#define BOOST_TEST_MODULE "miscellaneous"


#include "miscellaneous.hpp"

#include <boost/math/constants/constants.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( jacobiRotationMatrix ) {

    // A random Jacobi matrix is unitary
    BOOST_CHECK(GQCP::jacobiRotationMatrix(GQCP::JacobiRotationParameters(7, 4, 6.9921), 10).isUnitary());
    BOOST_CHECK(GQCP::jacobiRotationMatrix(GQCP::JacobiRotationParameters(9, 1, 78.00166), 22).isUnitary());

    // Let's see if we can construct the easiest Jacobi matrix, one with theta = pi/2 and dimension 2
    // cos(pi/2) = 0, sin(pi/2) = 1
    auto pi = boost::math::constants::half_pi<double>();
    Eigen::MatrixXd J = GQCP::jacobiRotationMatrix(GQCP::JacobiRotationParameters(1, 0, pi), 2);

    BOOST_CHECK(std::abs(J(0,0) - 0) < 1.0e-12);
    BOOST_CHECK(std::abs(J(0,1) - (-1)) < 1.0e-12);
    BOOST_CHECK(std::abs(J(1,0) - 1) < 1.0e-12);
    BOOST_CHECK(std::abs(J(0,0) - 0) < 1.0e-12);
}
