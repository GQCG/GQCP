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

#define BOOST_TEST_MODULE "SquareRankFourTensor"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/SquareRankFourTensor.hpp"


/**
 *  Check if the handling of the constructor argument behaves as expected.
 */
BOOST_AUTO_TEST_CASE(square_constructor) {

    GQCP::Tensor<double, 4> T1 {2, 2, 2, 2};
    BOOST_CHECK_NO_THROW(GQCP::SquareRankFourTensor<double> square_T1 {T1});

    GQCP::Tensor<double, 4> T2 {2, 1, 2, 2};
    BOOST_CHECK_THROW(GQCP::SquareRankFourTensor<double> square_T2 {T2}, std::invalid_argument);  // Not 'square'.
}


/**
 *  Check the error handling for `FromFile`.
 */
BOOST_AUTO_TEST_CASE(readArrayFromFile_tensor_throw) {

    const size_t dim = 7;

    // Check that there's an error when a wrong path is supplied.
    BOOST_CHECK_THROW(GQCP::SquareRankFourTensor<double>::FromFile("data/h2o_sto-3g_two_electron_horton.dat", dim), std::runtime_error);  // Should be 'datA'.


    // Check that there's no error when a correct path is supplied.
    BOOST_CHECK_NO_THROW(GQCP::SquareRankFourTensor<double>::FromFile("data/h2o_sto-3g_two_electron_horton.data", dim));


    // Check that there's an error when the tensor is incompatible with the given file.
    BOOST_CHECK_THROW(GQCP::SquareRankFourTensor<double>::FromFile("data/h2o_sto-3g_kinetic.data_horton", dim), std::runtime_error);  // Can't read in one-electron data in a tensor.
}


/**
 *  Check if `FromFile` method reads in a tensor correctly.
 */
BOOST_AUTO_TEST_CASE(readArrayFromFile_tensor_example) {

    // Test the read function on a small example mimicking the two-electron integrals.
    size_t dim = 6;
    GQCP::SquareRankFourTensor<double> T_ref {dim};
    T_ref.setZero();

    T_ref(0, 0, 0, 0) = 4.78506540471;
    T_ref(0, 0, 0, 1) = 0.741380351973;
    T_ref(0, 0, 0, 2) = 0.0;
    T_ref(0, 0, 0, 3) = 3.94054708595e-17;
    T_ref(0, 0, 0, 4) = 0.0;
    T_ref(0, 0, 0, 5) = 0.121785318177;
    T_ref(0, 0, 0, 6) = 0.121785318177;

    auto T = GQCP::SquareRankFourTensor<double>::FromFile("data/small_two_ints.data", dim);

    BOOST_CHECK(T.isApprox(T_ref, 1.0e-08));
}


/**
 *  Check if the pair-wise reshaping of a tensor into a matrix is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(pairWiseStrictReduced) {

    // Example 1.
    size_t dim1 = 2;
    GQCP::SquareRankFourTensor<double> T1 {dim1};
    for (size_t i = 0; i < dim1; i++) {
        for (size_t j = 0; j < dim1; j++) {
            for (size_t k = 0; k < dim1; k++) {
                for (size_t l = 0; l < dim1; l++) {
                    T1(i, j, k, l) = l + 2 * k + 4 * j + 8 * i;
                }
            }
        }
    }


    GQCP::SquareMatrix<double> M1_ref {1};  // 2*(2-1)/2 = 1
    M1_ref << 10;                           // By manual inspection we find that this is the only value that should appear in the matrix.

    BOOST_CHECK(M1_ref.isApprox(T1.pairWiseStrictReduced()));


    // Example 2.
    size_t dim2 = 3;
    GQCP::SquareRankFourTensor<double> T2 {dim2};
    for (size_t i = 0; i < dim2; i++) {
        for (size_t j = 0; j < dim2; j++) {
            for (size_t k = 0; k < dim2; k++) {
                for (size_t l = 0; l < dim2; l++) {
                    T2(i, j, k, l) = l + 3 * k + 9 * j + 27 * i;
                }
            }
        }
    }


    GQCP::SquareMatrix<double> M2_ref {3};  // 3*(3-1)/2 = 3
    // clang-format off
    M2_ref << 30, 33, 34,
              57, 60, 61,
              66, 69, 70;  // By manual inspection, these should be the elements of the reduced matrix.
    // clang-format on

    BOOST_CHECK(M2_ref.isApprox(T2.pairWiseStrictReduced()));
}
