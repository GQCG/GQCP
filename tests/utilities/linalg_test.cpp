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
#define BOOST_TEST_MODULE "linalg_test"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


#include "utilities/linalg.hpp"



BOOST_AUTO_TEST_CASE ( areEqual_throws ) {

    // Check for a throw if the dimensions aren't compatible.
    Eigen::Tensor<double, 4> M (2, 2, 2, 2);  // don't need to set them to zeros to check dimensions
    Eigen::Tensor<double, 4> T1 (2, 2, 3, 2);  // don't need to set them to zeros to check dimensions

    BOOST_CHECK_THROW(GQCP::areEqual(M, T1, 1.0e-6), std::invalid_argument);


    // Check for no throw if correct tensors are given.
    Eigen::Tensor<double, 4> T2 (2, 2, 2, 2);

    BOOST_CHECK_NO_THROW(GQCP::areEqual(M, T2, 1.0e-6));
}


BOOST_AUTO_TEST_CASE ( areEqual_example ) {

    Eigen::Tensor<double, 4> M (2, 2, 2, 2);
    Eigen::Tensor<double, 4> T (2, 2, 2, 2);


    // Fill the two compatible tensors with random data and check if they're not equal.
    M.setRandom();
    T.setRandom();

    BOOST_CHECK(!GQCP::areEqual(M, T, 1.0e-6));  // probability of random tensors being equal approaches zero


    // Fill the two compatible tensors with the same data and check if they're equal.
    M.setZero();
    T.setZero();
    M(0,1,0,0) = 0.5;
    T(0,1,0,0) = 0.5;

    BOOST_CHECK(GQCP::areEqual(M, T, 1.0e-6));
}


BOOST_AUTO_TEST_CASE ( areEqualEigenvectors ) {

    // Test areEqualEigenvectors with an example.
    Eigen::VectorXd a (3);
    Eigen::VectorXd b (3);
    Eigen::VectorXd c (3);
    Eigen::VectorXd d (3);

    a << 2, 3, 1;
    b << 2, 3, 1;
    c << -2, -3, -1;
    d << 2, 3, 0;


    BOOST_CHECK(GQCP::areEqualEigenvectors(a, b, 1.0e-6));
    BOOST_CHECK(GQCP::areEqualEigenvectors(a, c, 1.0e-6));
    BOOST_CHECK(GQCP::areEqualEigenvectors(b, c, 1.0e-6));

    BOOST_CHECK(!GQCP::areEqualEigenvectors(a, d, 1.0e-6));
    BOOST_CHECK(!GQCP::areEqualEigenvectors(c, d, 1.0e-6));
}


BOOST_AUTO_TEST_CASE ( areEqualSetsOfEigenvectors_throws ) {

    // Check for throws if the dimensions aren't compatible.
    Eigen::MatrixXd C1 (3, 3);
    Eigen::MatrixXd C2 (3, 2);

    BOOST_CHECK_THROW(GQCP::areEqualSetsOfEigenvectors(C1, C2, 1.0e-6), std::invalid_argument);


    // Check for no throw if the dimensions are compatible
    Eigen::MatrixXd C3 (3, 3);

    BOOST_CHECK_NO_THROW(GQCP::areEqualSetsOfEigenvectors(C1, C3, 1.0e-6));
}


BOOST_AUTO_TEST_CASE ( areEqualSetsOfEigenvectors_example ) {

    // Test areEqualSetsOfEigenvectors with an example
    Eigen::MatrixXd eigenvectors1 (2, 2);
    Eigen::MatrixXd eigenvectors2 (2, 2);
    Eigen::MatrixXd eigenvectors3 (2, 2);
    Eigen::MatrixXd eigenvectors4 (2, 2);

    eigenvectors1 << 0,  2,
                     1, -1;

    eigenvectors2 << 0,  2,
                     1, -1;

    eigenvectors3 << 0, -2,
                     1,  1;

    eigenvectors4 << 1,  2,
                     1, -1;


    BOOST_CHECK(GQCP::areEqualSetsOfEigenvectors(eigenvectors1, eigenvectors2, 1.0e-6));
    BOOST_CHECK(GQCP::areEqualSetsOfEigenvectors(eigenvectors1, eigenvectors3, 1.0e-6));

    BOOST_CHECK(!GQCP::areEqualSetsOfEigenvectors(eigenvectors1, eigenvectors4, 1.0e-6));
}


BOOST_AUTO_TEST_CASE ( strictLowerTriangle_matrix_invalid ) {

    Eigen::MatrixXd A_invalid (3, 4);
    BOOST_CHECK_THROW(GQCP::strictLowerTriangle(A_invalid), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( strictLowerTriangle_matrix ) {

    Eigen::MatrixXd A (3, 3);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    Eigen::VectorXd ref_strict_lower_triangle_A (3);
    ref_strict_lower_triangle_A << 4, 7, 8;

    BOOST_CHECK(ref_strict_lower_triangle_A.isApprox(GQCP::strictLowerTriangle(A)));


    Eigen::MatrixXd B (4, 4);
    B <<  1,  2,  3,  4,
          5,  6,  7,  8,
          9, 10, 11, 12,
         13, 14, 15, 16;

    Eigen::VectorXd ref_strict_lower_triangle_B (6);
    ref_strict_lower_triangle_B << 5, 9, 13, 10, 14, 15;

    BOOST_CHECK(ref_strict_lower_triangle_B.isApprox(GQCP::strictLowerTriangle(B)));
}


BOOST_AUTO_TEST_CASE ( fillStrictLowerTriangle_invalid ) {

    Eigen::VectorXd a_invalid (4);
    BOOST_CHECK_THROW(GQCP::fillStrictLowerTriangle(a_invalid), std::invalid_argument);  // 4 is not a triangular number
}


BOOST_AUTO_TEST_CASE ( fillStrictLowerTriangle ) {

    Eigen::VectorXd a (3);
    a << 1, 2, 3;

    Eigen::MatrixXd A_ref (3, 3);
    A_ref << 0, 0, 0,
             1, 0, 0,
             2, 3, 0;

    BOOST_CHECK(A_ref.isApprox(GQCP::fillStrictLowerTriangle(a), 1.0e-12));


    Eigen::VectorXd b (6);
    b << 1, 2, 3, 4, 5, 6;

    Eigen::MatrixXd B_ref (4, 4);
    B_ref << 0, 0, 0, 0,
             1, 0, 0, 0,
             2, 4, 0, 0,
             3, 5, 6, 0;

    BOOST_CHECK(B_ref.isApprox(GQCP::fillStrictLowerTriangle(b), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( toMatrix ) {

    // Create an example 2x2x2x2 tensor
    Eigen::Tensor<double, 4> T1 (2, 2, 2, 2);

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    T1(i,j,k,l) = l + 2*k + 4*j + 8*i;
                }
            }
        }
    }

    Eigen::MatrixXd M1_ref (4, 4);
    M1_ref <<  0,  2,  1,  3,
               8, 10,  9, 11,
               4,  6,  5,  7,
              12, 14, 13, 15;

    BOOST_CHECK(M1_ref.isApprox(GQCP::toMatrix(T1), 1.0e-12));


    // Create an example 3x3x3x3 tensor
    Eigen::Tensor<double, 4> T2 (3, 3, 3, 3);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
                for (size_t l = 0; l < 3; l++) {
                    T2(i,j,k,l) = l + 3*k + 9*j + 27*i;
                }
            }
        }
    }

    Eigen::MatrixXd M2_ref (9, 9);
    M2_ref <<  0,  3,  6,  1,  4,  7,  2,  5,  8,
              27, 30, 33, 28, 31, 34, 29, 32, 35,
              54, 57, 60, 55, 58, 61, 56, 59, 62,
               9, 12, 15, 10, 13, 16, 11, 14, 17,
              36, 39, 42, 37, 40, 43, 38, 41, 44,
              63, 66, 69, 64, 67, 70, 65, 68, 71,
              18, 21, 24, 19, 22, 25, 20, 23, 26,
              45, 48, 51, 46, 49, 52, 47, 50, 53,
              72, 75, 78, 73, 76, 79, 74, 77, 80;

    BOOST_CHECK(M2_ref.isApprox(GQCP::toMatrix(T2), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( strictLowerTriangle_tensor_invalid ) {

    Eigen::Tensor<double, 4> T_invalid (4, 3, 3, 3);
    BOOST_CHECK_THROW(GQCP::strictLowerTriangle(T_invalid), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( strictLowerTriangle_tensor ) {

    // Create an example tensor
    Eigen::Tensor<double, 4> T1 (2, 2, 2, 2);

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    T1(i,j,k,l) = l + 2*k + 4*j + 8*i;
                }
            }
        }
    }


    Eigen::MatrixXd M1_ref (1, 1);  // 2*(2-1)/2 = 1
    M1_ref << 10;  // by manual inspection we find that this is the only value that should appear in the matrix

    BOOST_CHECK(M1_ref.isApprox(GQCP::strictLowerTriangle(T1)));


    Eigen::Tensor<double, 4> T2 (3, 3, 3, 3);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
                for (size_t l = 0; l < 3; l++) {
                    T2(i,j,k,l) = l + 3*k + 9*j + 27*i;
                }
            }
        }
    }


    Eigen::MatrixXd M2_ref (3, 3);  // 3*(3-1)/2 = 3
    M2_ref << 30, 33, 34,
              57, 60, 61,
              66, 69, 70;  // by manual inspection, these should be the elements of the reduced matrix

    BOOST_CHECK(M2_ref.isApprox(GQCP::strictLowerTriangle(T2)));
}
