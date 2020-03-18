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
#define BOOST_TEST_MODULE "USQOneElectronOperator"

#include <boost/test/unit_test.hpp>
#include <math.h>

#include "Basis/ScalarBasis/CartesianGTO.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"

#include "Utilities/miscellaneous.hpp"


/**
 *  Check the construction of one-electron operators from matrices.
 */
BOOST_AUTO_TEST_CASE ( USQOneElectronOperator_constructor ) {

    // Check a correct constructor.
    const auto square_matrix = GQCP::SquareMatrix<double>::Zero(4, 4);
    const GQCP::ScalarUSQOneElectronOperator<double> O (square_matrix, square_matrix); 


    // Check a faulty constructor.
    // Either both the alpha and beta matrix dimensions are wrong, or only one of them is. All cases are checked.
    const GQCP::MatrixX<double> matrix = GQCP::MatrixX<double>::Zero(3, 4);
    BOOST_CHECK_THROW(GQCP::ScalarUSQOneElectronOperator<double> O_2 (matrix, matrix), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::ScalarUSQOneElectronOperator<double> O_3 (matrix, square_matrix), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::ScalarUSQOneElectronOperator<double> O_4 (square_matrix, matrix), std::invalid_argument);
}


/**
 *  Check if the zero constructor actually sets its parameters to zeros
 */
BOOST_AUTO_TEST_CASE ( USQOneElectronOperator_zero_constructor ) {

    const size_t dim = 2;
    const GQCP::ScalarUSQOneElectronOperator<double> one_op (2, 2);  // should initialize to zeros

    BOOST_CHECK_EQUAL(one_op.alphaDimension(), dim);
    BOOST_CHECK_EQUAL(one_op.betaDimension(), dim);
    BOOST_CHECK(one_op.alphaParameters().isZero(1.0e-08));
    BOOST_CHECK(one_op.betaParameters().isZero(1.0e-08));
}


/**
 *  Check if addition of operators works as expected
 */
BOOST_AUTO_TEST_CASE ( USQOneElectronOperator_addition ) {

    const size_t dim = 2;

    // Initialize two test matrices and convert them into operators
    GQCP::QCMatrix<double> M1 (dim);
    M1 << 1.0, 2.0,
          3.0, 4.0;
    const GQCP::ScalarUSQOneElectronOperator<double> op1 (M1, M1);

    GQCP::QCMatrix<double> M2 (dim);
    M2 << 5.0, 6.0,
          7.0, 8.0;
    const GQCP::ScalarUSQOneElectronOperator<double> op2 (M2, M2);


    // Initialize the reference and check the result
    GQCP::QCMatrix<double> M_sum_ref (dim);
    M_sum_ref <<  6.0,  8.0,
                 10.0, 12.0;
    
    const auto op_sum = op1 + op2;
    BOOST_CHECK(op_sum.alphaParameters().isApprox(M_sum_ref, 1.0e-08));
    BOOST_CHECK(op_sum.betaParameters().isApprox(M_sum_ref, 1.0e-08));
}


/**
 *  Check if difference of operators works as expected
 */
BOOST_AUTO_TEST_CASE ( USQOneElectronOperator_difference ) {

    const size_t dim = 2;

    // Initialize two test matrices and convert them into operators
    GQCP::QCMatrix<double> M1 (dim);
    M1 << 1.0, 2.0,
          3.0, 4.0;
    const GQCP::ScalarUSQOneElectronOperator<double> op1 (M1, M1);

    GQCP::QCMatrix<double> M2 (dim);
    M2 << 5.0, 6.0,
          7.0, 8.0;
    const GQCP::ScalarUSQOneElectronOperator<double> op2 (M2, M2);


    // Initialize the reference and check the result
    GQCP::QCMatrix<double> M_diff_ref (dim);
    M_diff_ref <<  4.0,  4.0,
                   4.0,  4.0;
    
    const auto op_diff = op2 - op1;
    BOOST_CHECK(op_diff.alphaParameters().isApprox(M_diff_ref, 1.0e-08));
    BOOST_CHECK(op_diff.betaParameters().isApprox(M_diff_ref, 1.0e-08));
}


/**
 *  Check if scalar product with the operators works as expected
 */
BOOST_AUTO_TEST_CASE ( USQOneElectronOperator_scalar_product ) {

    const size_t dim = 2;
    const double scalar = 2.0;

    // Initialize a test matrix and convert it into an operator
    GQCP::QCMatrix<double> M1 (dim);
    M1 << 1.0, 2.0,
          3.0, 4.0;
    const GQCP::ScalarUSQOneElectronOperator<double> op1 (M1, M1);

    // Initialize the reference and check the result
    GQCP::QCMatrix<double> M_ref (dim);
    M_ref <<  2.0,  4.0,
              6.0,  8.0;
    
    const auto op_result = scalar * op1;
    BOOST_CHECK(op_result.alphaParameters().isApprox(M_ref, 1.0e-08));
    BOOST_CHECK(op_result.betaParameters().isApprox(M_ref, 1.0e-08));
}


/**
 *  Check if negation of the operators works as expected
 */
BOOST_AUTO_TEST_CASE ( USQOneElectronOperator_negate ) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator
    GQCP::QCMatrix<double> M1 (dim);
    M1 << 1.0, 2.0,
          3.0, 4.0;
    const GQCP::ScalarUSQOneElectronOperator<double> op1 (M1, M1);

    // Initialize the reference and check the result
    GQCP::QCMatrix<double> M_ref (dim);
    M_ref <<  -1.0,  -2.0,
              -3.0,  -4.0;
    
    const auto op_result = -op1;
    BOOST_CHECK(op_result.alphaParameters().isApprox(M_ref, 1.0e-08));
    BOOST_CHECK(op_result.betaParameters().isApprox(M_ref, 1.0e-08));
}


/**
 *  Check if calculateExpectationValue throws when necessary
 */
BOOST_AUTO_TEST_CASE ( calculateExpectationValue_throw ) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator
    GQCP::QCMatrix<double> M1 (dim);
    M1 << 0.0, 0.0,
          0.0, 0.0;

    const GQCP::ScalarUSQOneElectronOperator<double> h (M1, M1);
    const GQCP::OneRDM<double> D_valid = GQCP::OneRDM<double>::Zero(2, 2);
    const GQCP::OneRDM<double> D_invalid = GQCP::OneRDM<double>::Zero(3, 3);

    BOOST_CHECK_THROW(h.calculateExpectationValue(D_invalid, D_invalid), std::invalid_argument);
    BOOST_CHECK_THROW(h.calculateExpectationValue(D_invalid, D_valid), std::invalid_argument);
    BOOST_CHECK_THROW(h.calculateExpectationValue(D_valid, D_invalid), std::invalid_argument);

    BOOST_CHECK_NO_THROW(h.calculateExpectationValue(D_valid, D_valid));
}


/**
 * Check whether or not calculateExpectationValue shows the correct behaviour
 */
BOOST_AUTO_TEST_CASE ( calculateExpectationValue_behaviour ) {
    
    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator
    GQCP::QCMatrix<double> M1 (dim);
    M1 << 1.0, 2.0,
          3.0, 4.0;
    const GQCP::ScalarUSQOneElectronOperator<double> op (M1, M1);

    // initialize an alpha and beta density matrix, each one is chosen as a hermitian matrix.
    GQCP::QCMatrix<double> d_alpha (dim);
    d_alpha << 0.0, 1.0,
               1.0, 0.0;

    GQCP::QCMatrix<double> d_beta (dim);
    d_beta << 1.0, 0.0,
              0.0, -1.0;
    
    // Initialize a reference value
    GQCP::QCMatrix<double> ref (1);
    ref << 2.0;

    const auto ex_value = op.calculateExpectationValue(d_alpha, d_beta);
    BOOST_CHECK(ex_value.isApprox(ref, 1.0e-08));
}


/**
 * Check whether or not the rotate with transformation matrix method works as expected
 */
BOOST_AUTO_TEST_CASE ( rotate_with_unitary_transformation_matrix ) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator
    GQCP::QCMatrix<double> M1 (dim);
    M1 << 1.0, 2.0,
          3.0, 4.0;
    GQCP::ScalarUSQOneElectronOperator<double> op (M1, M1);

    // Initialize a unitary transformation matrix
    GQCP::TransformationMatrix<double> U (dim);
    U << 1.0, 0.0,
         0.0, 1.0;
    
    op.rotate(U);
    BOOST_CHECK(op.alphaParameters().isApprox(M1, 1.0e-08));
    BOOST_CHECK(op.betaParameters().isApprox(M1, 1.0e-08));

}


/**
 * Check whether or not the rotate with transformation matrix method works as expected
 */
BOOST_AUTO_TEST_CASE ( transform_with_transformation_matrix ) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator
    GQCP::QCMatrix<double> M1 (dim);
    M1 << 1.0, 2.0,
          3.0, 4.0;
    GQCP::ScalarUSQOneElectronOperator<double> op (M1, M1);

    // Initialize a transformation matrix
    GQCP::TransformationMatrix<double> T (dim);
    T << 1.0, 2.0,
         3.0, 4.0;

    // Initialize a reference matrix
    GQCP::QCMatrix<double> ref (dim);
    ref << 7.0, 10,
           15.0, 22.0;
    
    op.transform(T);
    BOOST_CHECK(op.alphaParameters().isApprox(ref, 1.0e-08));
    BOOST_CHECK(op.betaParameters().isApprox(ref, 1.0e-08));
}


/**
 * Check whether or not the rotate with transformation matrix method works as expected
 */
BOOST_AUTO_TEST_CASE ( transform_with_jacobi_matrix ) {

    const size_t dim = 3;

    // Initialize a test matrix and convert it into an operator
    GQCP::QCMatrix<double> M1 (dim);
    M1 << 1.0, 7.0, 3.0,
          7.0, 4.0, -5.0,
          3.0, -5.0, 6.0;
    GQCP::ScalarUSQOneElectronOperator<double> op (M1, M1);

    // Initialize a transformation matrix
    GQCP::JacobiRotationParameters J (1, 2, M_PI_2);

    // Initialize a reference matrix
    GQCP::QCMatrix<double> ref (dim);
    ref <<  -4.65891, 0.00000, 5.47639,
            0.00000, 9.65891, -2.00230,
            5.47639, -2.00230, 6.00000;

    op.rotate(J);
    BOOST_CHECK(op.alphaParameters().isApprox(ref, 1.0e-08));
    BOOST_CHECK(op.betaParameters().isApprox(ref, 1.0e-08));
}
