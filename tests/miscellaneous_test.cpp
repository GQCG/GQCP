// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
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


BOOST_AUTO_TEST_CASE ( minors ) {

    Eigen::MatrixXd A (3, 4);
    A << 1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12;


    Eigen::MatrixXd A_00 (2, 3);
    A_00 <<  6,  7,  8,
            10, 11, 12;
    BOOST_CHECK(A_00.isApprox(GQCP::minor(A, 0, 0)));

    Eigen::MatrixXd A_21 (2, 3);
    A_21 << 1, 3, 4,
            5, 7, 8;
    BOOST_CHECK(A_21.isApprox(GQCP::minor(A, 2, 1)));
}


BOOST_AUTO_TEST_CASE ( permanent_throws ) {

    Eigen::MatrixXd A (3, 4);
    Eigen::MatrixXd B (3, 2);

    BOOST_CHECK_THROW(GQCP::permanent(A), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::permanent(B), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( permanent ) {

    Eigen::MatrixXd A (2, 2);
    A << 2, 3,
         9, 1;
    BOOST_CHECK(std::abs(GQCP::permanent(A) - 29.0) < 1.0e-12);


    Eigen::MatrixXd B (3, 3);
    B << 1,  2, -3,
         4, -5,  6,
         7, -8,  9;
    BOOST_CHECK(std::abs(GQCP::permanent(B) - 264.0) < 1.0e-12);

}
