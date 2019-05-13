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
#define BOOST_TEST_MODULE "SquareRankFourTensor"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "math/SquareRankFourTensor.hpp"


BOOST_AUTO_TEST_CASE ( square_constructor ) {

    GQCP::Tensor<double, 4> T1 (2, 2, 2, 2);
    T1.setZero();
    BOOST_CHECK_NO_THROW(GQCP::SquareRankFourTensor<double> square_T1 (T1));

    GQCP::Tensor<double, 4> T2 (2, 1, 2, 2);
    BOOST_CHECK_THROW(GQCP::SquareRankFourTensor<double> square_T2 (T2), std::invalid_argument);  // not square
}


BOOST_AUTO_TEST_CASE ( constructor_assignment ) {

    // A small check to see if the interface of the constructor and assignment operator works as expected

    GQCP::SquareRankFourTensor<double> A (2);
    GQCP::SquareRankFourTensor<double> B (2);

    GQCP::SquareRankFourTensor<double> T1 (A.Eigen() + B.Eigen());
    GQCP::SquareRankFourTensor<double> T2 = 2 * B.Eigen();

    GQCP::SquareRankFourTensor<double> T3 (T1.Eigen() + T2.Eigen());
    GQCP::SquareRankFourTensor<double> T4 = (T1.Eigen() + T2.Eigen());
    GQCP::SquareRankFourTensor<double> T5 = 3 * T2.Eigen();
}


BOOST_AUTO_TEST_CASE ( readArrayFromFile_tensor_throw ) {

    size_t dim = 7;

    // Check that there's an error when a wrong path is supplied
    BOOST_CHECK_THROW(GQCP::SquareRankFourTensor<double>::FromFile("data/h2o_sto-3g_two_electron_horton.dat", dim), std::runtime_error);


    // Check that there's no error when a correct path is supplied
    BOOST_CHECK_NO_THROW(GQCP::SquareRankFourTensor<double>::FromFile("data/h2o_sto-3g_two_electron_horton.data", dim));


    // Check that there's an error when the tensor is incompatible with the given file
    BOOST_CHECK_THROW(GQCP::SquareRankFourTensor<double>::FromFile("data/h2o_sto-3g_kinetic.data_horton", dim), std::runtime_error);  // can't read in one-electron data in a tensor
}


BOOST_AUTO_TEST_CASE ( readArrayFromFile_tensor_example ) {

    // Test the read function on a small example mimicking the two-electron integrals
    size_t dim = 6;
    GQCP::SquareRankFourTensor<double> T_ref (dim);
    T_ref.setZero();

    T_ref(0,0,0,0) = 4.78506540471;
    T_ref(0,0,0,1) = 0.741380351973;
    T_ref(0,0,0,2) = 0.0;
    T_ref(0,0,0,3) = 3.94054708595e-17;
    T_ref(0,0,0,4) = 0.0;
    T_ref(0,0,0,5) = 0.121785318177;
    T_ref(0,0,0,6) = 0.121785318177;

    auto T = GQCP::SquareRankFourTensor<double>::FromFile("data/small_two_ints.data", dim);

    BOOST_CHECK(T.isApprox(T_ref, 1.0e-08));
}


BOOST_AUTO_TEST_CASE ( pairWiseStrictReduce ) {

    // Example 1
    size_t dim1 = 2;
    GQCP::SquareRankFourTensor<double> T1 (dim1);
    for (size_t i = 0; i < dim1; i++) {
        for (size_t j = 0; j < dim1; j++) {
            for (size_t k = 0; k < dim1; k++) {
                for (size_t l = 0; l < dim1; l++) {
                    T1(i,j,k,l) = l + 2*k + 4*j + 8*i;
                }
            }
        }
    }


    GQCP::SquareMatrix<double> M1_ref (1);  // 2*(2-1)/2 = 1
    M1_ref << 10;  // by manual inspection we find that this is the only value that should appear in the matrix

    BOOST_CHECK(M1_ref.isApprox(T1.pairWiseStrictReduce()));


    // Example 2
    size_t dim2 = 3;
    GQCP::SquareRankFourTensor<double> T2 (dim2);
    for (size_t i = 0; i < dim2; i++) {
        for (size_t j = 0; j < dim2; j++) {
            for (size_t k = 0; k < dim2; k++) {
                for (size_t l = 0; l < dim2; l++) {
                    T2(i,j,k,l) = l + 3*k + 9*j + 27*i;
                }
            }
        }
    }


    GQCP::SquareMatrix<double> M2_ref (3);  // 3*(3-1)/2 = 3
    M2_ref << 30, 33, 34,
              57, 60, 61,
              66, 69, 70;  // by manual inspection, these should be the elements of the reduced matrix

    BOOST_CHECK(M2_ref.isApprox(T2.pairWiseStrictReduce()));
}


BOOST_AUTO_TEST_CASE ( SquareRankFourTensor_matrix_contraction_trivial ) {

    // toy rank-4 tensor
    GQCP::SquareRankFourTensor<double> T (2);
    for (long i = 0; i < 2; i++) {
        for (long j = 0; j < 2; j++) {
            for (long k = 0; k < 2; k++) {
                for (long l = 0; l < 2; l++) {
                    T(i, j, k, l) = l + 2 * k + 4 * j + 8 * i;
                }
            }
        }
    }

    // shuffled tensor
    std::array<int, 4> shuffle {3,2,1,0};
    GQCP::SquareRankFourTensor<double> T2  =  GQCP::SquareRankFourTensor<double>(T.Eigen().shuffle(shuffle));

    // toy matrix
    GQCP::SquareMatrix<double> A (2);
    for (long i = 0; i < 2; i++) {
        for (long j = 0; j < 2; j++) {
            A(i, j) = 1 + j + 2 * i;
        }
    }

    // test rank-4 tensor, data from https://stackoverflow.com/questions/47556726
    GQCP::SquareRankFourTensor<double> T_test = GQCP::SquareRankFourTensor<double>::FromFile("data/tensor_contraction_test.data", 2);

    T.contract(A, 0);

    BOOST_CHECK(T_test.isApprox(T_test, 1.0e-12));

    // If we contract the with the last index axis of the shuffled tensor, then shuffle back the results should be the same.
    T2.contract(A, 3);
    GQCP::SquareRankFourTensor<double> T3 = GQCP::SquareRankFourTensor<double>(T2.Eigen().shuffle(shuffle));

    BOOST_CHECK(T3.isApprox(T_test, 1.0e-12));

    BOOST_CHECK_THROW(T3.contract(A, 5), std::invalid_argument);
}
