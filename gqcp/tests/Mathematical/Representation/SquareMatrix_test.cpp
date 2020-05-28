// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE "SquareMatrix"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/SquareMatrix.hpp"


BOOST_AUTO_TEST_CASE(constructor) {

    auto M1 = GQCP::SquareMatrix<double>::Zero(2, 2);
    BOOST_CHECK_NO_THROW(GQCP::SquareMatrix<double> square_M1 {M1});

    auto M2 = GQCP::SquareMatrix<double>::Zero(2, 1);
    BOOST_CHECK_THROW(GQCP::SquareMatrix<double> square_M2 {M2}, std::invalid_argument);
}


BOOST_AUTO_TEST_CASE(constructor_assignment) {

    // A small check to see if the interface of the constructor and assignment operator works as expected

    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Random(3, 3);
    GQCP::SquareMatrix<double> B = GQCP::SquareMatrix<double>::Random(3, 3);

    GQCP::SquareMatrix<double> M1 {A * B};
    GQCP::SquareMatrix<double> M2 = A + B;
    GQCP::SquareMatrix<double> M3 = 2 * A;

    GQCP::SquareMatrix<double> M4 {M1 * M2};
    GQCP::SquareMatrix<double> M5 = M2 + M3;
    GQCP::SquareMatrix<double> M6 = 4 * M2;
}


BOOST_AUTO_TEST_CASE(FromStrictTriangle_invalid) {

    GQCP::VectorX<double> a_invalid {4};
    BOOST_CHECK_THROW(GQCP::SquareMatrix<double>::FromStrictTriangle(a_invalid), std::invalid_argument);  // 4 is not a triangular number
}


BOOST_AUTO_TEST_CASE(FromStrictTriangle) {

    GQCP::VectorX<double> a {3};
    a << 1, 2, 3;

    GQCP::SquareMatrix<double> A_ref {3};
    // clang-format off
    A_ref << 0, 0, 0,
             1, 0, 0,
             2, 3, 0;
    // clang-format on

    BOOST_CHECK(A_ref.isApprox(GQCP::SquareMatrix<double>::FromStrictTriangle(a), 1.0e-12));


    GQCP::VectorX<double> b {6};
    b << 1, 2, 3, 4, 5, 6;

    GQCP::SquareMatrix<double> B_ref {4};
    // clang-format off
    B_ref << 0, 0, 0, 0,
             1, 0, 0, 0,
             2, 4, 0, 0,
             3, 5, 6, 0;
    // clang-format on

    BOOST_CHECK(B_ref.isApprox(GQCP::SquareMatrix<double>::FromStrictTriangle(b), 1.0e-12));
}


BOOST_AUTO_TEST_CASE(FromTriangle) {

    GQCP::VectorX<double> upper_triangle {6};
    upper_triangle << 1, 2, 3, 4, 5, 6;

    GQCP::SquareMatrix<double> H_ref {3};
    // clang-format off
    H_ref << 1, 2, 3,
             2, 4, 5,
             3, 5, 6;
    // clang-format on

    BOOST_CHECK(H_ref.isApprox(GQCP::SquareMatrix<double>::SymmetricFromUpperTriangle(upper_triangle)));
}


BOOST_AUTO_TEST_CASE(pairWiseStrictReduced) {

    GQCP::SquareMatrix<double> A {3};
    // clang-format off
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    // clang-format on

    GQCP::VectorX<double> ref_strict_lower_triangle_A {3};
    ref_strict_lower_triangle_A << 4, 7, 8;

    BOOST_CHECK(ref_strict_lower_triangle_A.isApprox(A.pairWiseStrictReduced()));


    GQCP::SquareMatrix<double> B {4};
    // clang-format off
    B <<  1,  2,  3,  4,
          5,  6,  7,  8,
          9, 10, 11, 12,
         13, 14, 15, 16;
    // clang-format on

    GQCP::VectorX<double> ref_strict_lower_triangle_B {6};
    ref_strict_lower_triangle_B << 5, 9, 13, 10, 14, 15;

    BOOST_CHECK(ref_strict_lower_triangle_B.isApprox(B.pairWiseStrictReduced()));
}


BOOST_AUTO_TEST_CASE(calculatePermanentCombinatorial) {

    GQCP::SquareMatrix<double> A {2};
    // clang-format off
    A << 2, 3,
         9, 1;
    // clang-format on
    BOOST_CHECK(std::abs(A.calculatePermanentCombinatorial() - 29.0) < 1.0e-12);


    GQCP::SquareMatrix<double> B {3};
    // clang-format off
    B << 1,  2, -3,
         4, -5,  6,
         7, -8,  9;
    // clang-format on
    BOOST_CHECK(std::abs(B.calculatePermanentCombinatorial() - 264.0) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE(calculatePermanentRyser) {

    GQCP::SquareMatrix<double> A {2};
    // clang-format off
    A << 2, 3,
         9, 1;
    // clang-format on
    BOOST_CHECK(std::abs(A.calculatePermanentRyser() - 29.0) < 1.0e-12);


    GQCP::SquareMatrix<double> B {3};
    // clang-format off
    B << 1,  2, -3,
         4, -5,  6,
         7, -8,  9;
    // clang-format on
    BOOST_CHECK(std::abs(B.calculatePermanentRyser() - 264.0) < 1.0e-12);


    GQCP::SquareMatrix<double> C = GQCP::SquareMatrix<double>::Random(5, 5);
    BOOST_CHECK(std::abs(C.calculatePermanentCombinatorial() - C.calculatePermanentRyser()) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE(noPivotLUDecompose) {

    // Reference data from https://stackoverflow.com/questions/41150997/perform-lu-decomposition-without-pivoting-in-matlab
    GQCP::SquareMatrix<double> A {3};
    GQCP::SquareMatrix<double> L_ref {3};
    GQCP::SquareMatrix<double> U_ref {3};
    // clang-format off
    A << 7, 6, 10,
         3, 8, 7,
         3, 5, 5;
    
    L_ref <<  1.0000, 0,      0,
              0.4286, 1,      0,
              0.4286, 0.4474, 1;

    U_ref <<  7, 6.0000,     10,
              0, 5.4286, 2.7143,
              0, 0,     -0.5000;
    // clang-format on

    const auto& LU = A.noPivotLUDecompose();

    const auto& L = LU[0];
    const auto& U = LU[1];

    BOOST_CHECK(U.isApprox(U_ref, 1.0e-4));
    BOOST_CHECK(L.isApprox(L_ref, 1.0e-4));


    // Reference data from https://www.geeksforgeeks.org/doolittle-algorithm-lu-decomposition/
    GQCP::SquareMatrix<double> A2 {3};
    GQCP::SquareMatrix<double> L_ref2 {3};
    GQCP::SquareMatrix<double> U_ref2 {3};
    // clang-format off
    A2 << 2, -1, -2,
         -4,  6,  3,
         -4, -2,  8;
    
    L_ref2 <<  1,  0,  0,
              -2,  1,  0,
              -2, -1,  1;

    U_ref2 << 2, -1, -2,
              0,  4, -1,
              0,  0,  3;
    // clang-format on


    const auto& LU2 = A2.noPivotLUDecompose();

    const auto& L2 = LU2[0];
    const auto& U2 = LU2[1];

    BOOST_CHECK(U2.isApprox(U_ref2, 1.0e-5));
    BOOST_CHECK(L2.isApprox(L_ref2, 1.0e-5));
}
