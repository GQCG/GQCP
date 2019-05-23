// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#define BOOST_TEST_MODULE "SquareMatrix"

#include <boost/test/unit_test.hpp>

#include "math/SquareMatrix.hpp"

#include <boost/math/constants/constants.hpp>


BOOST_AUTO_TEST_CASE ( constructor ) {

    auto M1 = GQCP::SquareMatrix<double>::Zero(2, 2);
    BOOST_CHECK_NO_THROW(GQCP::SquareMatrix<double> square_M1 (M1));

    auto M2 = GQCP::SquareMatrix<double>::Zero(2, 1);
    BOOST_CHECK_THROW(GQCP::SquareMatrix<double> square_M2 (M2), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( constructor_assignment ) {

    // A small check to see if the interface of the constructor and assignment operator works as expected

    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Random(3, 3);
    GQCP::SquareMatrix<double> B = GQCP::SquareMatrix<double>::Random(3, 3);

    GQCP::SquareMatrix<double> M1 (A * B);
    GQCP::SquareMatrix<double> M2 = A + B;
    GQCP::SquareMatrix<double> M3 = 2*A;

    GQCP::SquareMatrix<double> M4 (M1 * M2);
    GQCP::SquareMatrix<double> M5 = M2 + M3;
    GQCP::SquareMatrix<double> M6 = 4*M2;
}


BOOST_AUTO_TEST_CASE ( FromStrictTriangle_invalid ) {

    GQCP::VectorX<double> a_invalid (4);
    BOOST_CHECK_THROW(GQCP::SquareMatrix<double>::FromStrictTriangle(a_invalid), std::invalid_argument);  // 4 is not a triangular number
}


BOOST_AUTO_TEST_CASE ( FromStrictTriangle ) {

    GQCP::VectorX<double> a (3);
    a << 1, 2, 3;

    GQCP::SquareMatrix<double> A_ref (3);
    A_ref << 0, 0, 0,
             1, 0, 0,
             2, 3, 0;

    BOOST_CHECK(A_ref.isApprox(GQCP::SquareMatrix<double>::FromStrictTriangle(a), 1.0e-12));


    GQCP::VectorX<double> b (6);
    b << 1, 2, 3, 4, 5, 6;

    GQCP::SquareMatrix<double> B_ref (4);
    B_ref << 0, 0, 0, 0,
             1, 0, 0, 0,
             2, 4, 0, 0,
             3, 5, 6, 0;

    BOOST_CHECK(B_ref.isApprox(GQCP::SquareMatrix<double>::FromStrictTriangle(b), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( FromTriangle ) {

    GQCP::VectorX<double> upper_triangle (6);
    upper_triangle << 1, 2, 3, 4, 5, 6;

    GQCP::SquareMatrix<double> H_ref (3);
    H_ref << 1, 2, 3,
             2, 4, 5,
             3, 5, 6;

    BOOST_CHECK(H_ref.isApprox(GQCP::SquareMatrix<double>::FullFromTriangle(upper_triangle)));
}

//TODO: enable back after linking to GQCP 
//BOOST_AUTO_TEST_CASE ( FromJacobi ) {
//
//    // A random Jacobi matrix is unitary
//    BOOST_CHECK(GQCP::SquareMatrix<double>::FromJacobi(GQCP::JacobiRotationParameters(7, 4, 6.9921), 10).isUnitary());
//    BOOST_CHECK(GQCP::SquareMatrix<double>::FromJacobi(GQCP::JacobiRotationParameters(9, 1, 78.00166), 22).isUnitary());
//
//    // Let's see if we can construct the easiest Jacobi matrix, one with theta = pi/2 and dimension 2
//    // cos(pi/2) = 0, sin(pi/2) = 1
//    auto pi = boost::math::constants::half_pi<double>();
//    auto J = GQCP::SquareMatrix<double>::FromJacobi(GQCP::JacobiRotationParameters(1, 0, pi), 2);
//
//    BOOST_CHECK(std::abs(J(0,0) - 0) < 1.0e-12);
//    BOOST_CHECK(std::abs(J(0,1) - (-1)) < 1.0e-12);
//    BOOST_CHECK(std::abs(J(1,0) - 1) < 1.0e-12);
//    BOOST_CHECK(std::abs(J(0,0) - 0) < 1.0e-12);
//}


BOOST_AUTO_TEST_CASE ( strictLowerTriangle ) {

    GQCP::SquareMatrix<double> A (3);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    GQCP::VectorX<double> ref_strict_lower_triangle_A (3);
    ref_strict_lower_triangle_A << 4, 7, 8;

    BOOST_CHECK(ref_strict_lower_triangle_A.isApprox(A.strictLowerTriangle()));


    GQCP::SquareMatrix<double> B (4);
    B <<  1,  2,  3,  4,
          5,  6,  7,  8,
          9, 10, 11, 12,
         13, 14, 15, 16;

    GQCP::VectorX<double> ref_strict_lower_triangle_B (6);
    ref_strict_lower_triangle_B << 5, 9, 13, 10, 14, 15;

    BOOST_CHECK(ref_strict_lower_triangle_B.isApprox(B.strictLowerTriangle()));
}


BOOST_AUTO_TEST_CASE ( permanent_combinatorial ) {

    GQCP::SquareMatrix<double> A (2);
    A << 2, 3,
         9, 1;
    BOOST_CHECK(std::abs(A.permanent_combinatorial() - 29.0) < 1.0e-12);


    GQCP::SquareMatrix<double> B (3);
    B << 1,  2, -3,
         4, -5,  6,
         7, -8,  9;
    BOOST_CHECK(std::abs(B.permanent_combinatorial() - 264.0) < 1.0e-12);

}


BOOST_AUTO_TEST_CASE ( permanent_ryser ) {

    GQCP::SquareMatrix<double> A (2);
    A << 2, 3,
         9, 1;
    BOOST_CHECK(std::abs(A.permanent_ryser() - 29.0) < 1.0e-12);


    GQCP::SquareMatrix<double> B (3);
    B << 1,  2, -3,
         4, -5,  6,
         7, -8,  9;
    BOOST_CHECK(std::abs(B.permanent_ryser() - 264.0) < 1.0e-12);


    GQCP::SquareMatrix<double> C = GQCP::SquareMatrix<double>::Random(5, 5);
    BOOST_CHECK(std::abs(C.permanent_combinatorial() - C.permanent_ryser()) < 1.0e-12);
}
