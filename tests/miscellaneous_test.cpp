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

    BOOST_CHECK_THROW(GQCP::permanent_combinatorial(A), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::permanent_combinatorial(B), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( permanent ) {

    Eigen::MatrixXd A (2, 2);
    A << 2, 3,
         9, 1;
    BOOST_CHECK(std::abs(GQCP::permanent_combinatorial(A) - 29.0) < 1.0e-12);


    Eigen::MatrixXd B (3, 3);
    B << 1,  2, -3,
         4, -5,  6,
         7, -8,  9;
    BOOST_CHECK(std::abs(GQCP::permanent_combinatorial(B) - 264.0) < 1.0e-12);

}


BOOST_AUTO_TEST_CASE ( gray_code ) {

    BOOST_CHECK(GQCP::gray_code(0) == 0);    // "0000" (0)
    BOOST_CHECK(GQCP::gray_code(1) == 1);    // "0001" (1)
    BOOST_CHECK(GQCP::gray_code(2) == 3);    // "0011" (3)
    BOOST_CHECK(GQCP::gray_code(3) == 2);    // "0010" (2)
    BOOST_CHECK(GQCP::gray_code(4) == 6);    // "0110" (6)
    BOOST_CHECK(GQCP::gray_code(5) == 7);    // "0111" (7)
    BOOST_CHECK(GQCP::gray_code(6) == 5);    // "0101" (5)
    BOOST_CHECK(GQCP::gray_code(7) == 4);    // "0100" (4)
    BOOST_CHECK(GQCP::gray_code(8) == 12);   // "1100" (12)
    BOOST_CHECK(GQCP::gray_code(9) == 13);   // "1101" (13)
    BOOST_CHECK(GQCP::gray_code(10) == 15);  // "1111" (15)
    BOOST_CHECK(GQCP::gray_code(11) == 14);  // "1110" (14)
    BOOST_CHECK(GQCP::gray_code(12) == 10);  // "1010" (10)
    BOOST_CHECK(GQCP::gray_code(13) == 11);  // "1011" (11)
    BOOST_CHECK(GQCP::gray_code(14) == 9);   // "1001" (9)
    BOOST_CHECK(GQCP::gray_code(15) == 8);   // "1000" (8)
}


BOOST_AUTO_TEST_CASE ( permanent_ryser_throws ) {

    Eigen::MatrixXd A (3, 4);
    Eigen::MatrixXd B (3, 2);

    BOOST_CHECK_THROW(GQCP::permanent_ryser(A), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::permanent_ryser(B), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( permanent_ryser ) {

    Eigen::MatrixXd A (2, 2);
    A << 2, 3,
         9, 1;
    BOOST_CHECK(std::abs(GQCP::permanent_ryser(A) - 29.0) < 1.0e-12);


    Eigen::MatrixXd B (3, 3);
    B << 1,  2, -3,
         4, -5,  6,
         7, -8,  9;
    BOOST_CHECK(std::abs(GQCP::permanent_ryser(B) - 264.0) < 1.0e-12);


    Eigen::MatrixXd C = Eigen::MatrixXd::Random(5, 5);
    BOOST_CHECK(std::abs(GQCP::permanent_combinatorial(C) - GQCP::permanent_ryser(C)) < 1.0e-12);
}
