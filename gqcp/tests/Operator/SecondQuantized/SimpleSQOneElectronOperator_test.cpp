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

#define BOOST_TEST_MODULE "SimpleSQOneElectronOperator"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"


/**
 *  Check the construction of one-electron operators from matrices.
 * 
 *  Since RSQOneElectronOperator derives from SimpleSQOneElectronOperator, we test the base functionality using a derived class.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_constructor) {

    // Set up some test matrices to be used as one-electron operator integrals.
    const auto matrix34 = GQCP::MatrixX<double>::Zero(3, 4);

    const auto square_matrix33 = GQCP::SquareMatrix<double>::Zero(3, 3);
    const auto square_matrix44 = GQCP::SquareMatrix<double>::Zero(4, 4);


    // Check a correct constructor: a square matrix may represent the integrals of a one-electron operator.
    BOOST_CHECK_NO_THROW(GQCP::ScalarRSQOneElectronOperator<double> operator_square {square_matrix33});

    // Check a throwing constructor: rectangular matrices cannot represent the integrals of a one-electron operator.
    BOOST_CHECK_THROW(GQCP::ScalarRSQOneElectronOperator<double> operator_rectangular {matrix34}, std::invalid_argument);

    // Check a throwing constructor: the dimensions of the square matrices must be equal.
    BOOST_CHECK_THROW(GQCP::VectorRSQOneElectronOperator<double> operator_rectangular({square_matrix33, square_matrix33, square_matrix44}), std::invalid_argument);
}


/**
 *  Check if the zero constructor actually sets its parameters to zeros.
 * 
 *  Since RSQOneElectronOperator derives from SimpleSQOneElectronOperator, we test the base functionality using a derived class.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_zero_constructor) {

    const size_t dim = 2;

    // Check a zero constructor for scalar operators.
    const GQCP::ScalarRSQOneElectronOperator<double> zero_scalar_operator {2};

    BOOST_CHECK_EQUAL(zero_scalar_operator.numberOfOrbitals(), dim);
    BOOST_CHECK(zero_scalar_operator.parameters().isZero(1.0e-08));


    // Check a zero constructor for vector operators.
    const GQCP::VectorRSQOneElectronOperator<double> zero_vector_operator {2};

    BOOST_CHECK_EQUAL(zero_vector_operator.numberOfOrbitals(), dim);
    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(zero_vector_operator.parameters(i).isZero(1.0e-08));
    }
}


/**
 *  Check if we can access a single component of a SimpleSQOneElectronOperator easily using the operator().
 */
BOOST_AUTO_TEST_CASE(operator_call) {

    const size_t dim = 4;

    // Set up three (random) matrix representations and construct the corresponding vector operator.
    const GQCP::QCMatrix<double> M1 = GQCP::QCMatrix<double>::Random(dim, dim);
    const GQCP::QCMatrix<double> M2 = GQCP::QCMatrix<double>::Random(dim, dim);
    const GQCP::QCMatrix<double> M3 = GQCP::QCMatrix<double>::Random(dim, dim);

    const GQCP::VectorRSQOneElectronOperator<double> vector_operator {{M1, M2, M3}};

    // Check if we can access a single component of this operator.
    BOOST_CHECK(vector_operator(0).parameters().isApprox(M1, 1.0e-12));
    BOOST_CHECK(vector_operator(1).parameters().isApprox(M2, 1.0e-12));
    BOOST_CHECK(vector_operator(2).parameters().isApprox(M3, 1.0e-12));
}


/**
 *  Check if calculateExpectationValue throws when necessary.
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_throw) {

    // Initialize a test operator and density matrices.
    const GQCP::ScalarRSQOneElectronOperator<double> f_operator {2};
    const GQCP::OneDM<double> D_valid = GQCP::OneDM<double>::Zero(2, 2);
    const GQCP::OneDM<double> D_invalid = GQCP::OneDM<double>::Zero(3, 3);

    BOOST_CHECK_NO_THROW(f_operator.calculateExpectationValue(D_valid));
    BOOST_CHECK_THROW(f_operator.calculateExpectationValue(D_invalid), std::invalid_argument);
}


/**
 *  Check if calculateExpectationValue shows the correct behaviour.
 * 
 *  We're testing this behavior for 'restricted' operators, but since the function is a method on its base class `SimpleSQOneElectronOperator`, this test also checks the validity on 'general' operators.
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_behaviour) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator.
    GQCP::QCMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on
    const GQCP::ScalarRSQOneElectronOperator<double> op {M1};

    // Initialize a test density matrix.
    GQCP::QCMatrix<double> D {dim};
    // clang-format off
    D << 0.0, 1.0,
         1.0, 0.0;
    // clang-format on

    // Initialize a reference value and check the result.
    const double reference_expectation_value = 5.0;

    const auto expectation_value = op.calculateExpectationValue(D)();  // extract the 'scalar' from a one-dimensional array
    BOOST_CHECK(std::abs(expectation_value - reference_expectation_value) < 1.0e-08);
}


/**
 *  Check if the addition of operators works as expected.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_addition) {

    const size_t dim = 2;

    // Initialize two test matrices and convert them into operators.
    GQCP::QCMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    const GQCP::ScalarRSQOneElectronOperator<double> op1 {M1};
    // clang-format on

    GQCP::QCMatrix<double> M2 {dim};
    // clang-format off
    M2 << 5.0, 6.0,
          7.0, 8.0;
    // clang-format on
    const GQCP::ScalarRSQOneElectronOperator<double> op2 {M2};


    // Initialize the reference and check the result of the operator addition.
    GQCP::QCMatrix<double> M_sum_ref {dim};
    // clang-format off
    M_sum_ref <<  6.0,  8.0,
                 10.0, 12.0;
    // clang-format on

    const auto op_sum = op1 + op2;
    BOOST_CHECK(op_sum.parameters().isApprox(M_sum_ref, 1.0e-08));
}


/**
 *  Check if scalar multiplication of a one-electron operator works as expected.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_scalar_product) {

    const size_t dim = 2;
    const double scalar = 2.0;

    // Initialize a test matrix and convert it into an operator.
    GQCP::QCMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on
    const GQCP::ScalarRSQOneElectronOperator<double> op1 {M1};

    // Initialize the reference and check the result of the scalar multiplication.
    GQCP::QCMatrix<double> M_ref {dim};
    // clang-format off
    M_ref <<  2.0,  4.0,
              6.0,  8.0;
    // clang-format on

    const auto op_result1 = scalar * op1;
    BOOST_CHECK(op_result1.parameters().isApprox(M_ref, 1.0e-08));

    const auto op_result2 = op1 * scalar;
    BOOST_CHECK(op_result2.parameters().isApprox(M_ref, 1.0e-08));
}


/**
 *  Check if the negation of a one-electron operator works as expected.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_negation) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator.
    GQCP::QCMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on
    const GQCP::ScalarRSQOneElectronOperator<double> op1 {M1};

    // Initialize the reference and check the result of the negation.
    GQCP::QCMatrix<double> M_ref {dim};
    // clang-format off
    M_ref <<  -1.0,  -2.0,
              -3.0,  -4.0;
    // clang-format on

    const auto op_result = -op1;
    BOOST_CHECK(op_result.parameters().isApprox(M_ref, 1.0e-08));
}


/**
 *  Check if dot product multiplication of a SimpleSQOneElectronOperator with a Vector with the same components is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(dot) {

    // Set up a toy RSQOneElectronOperator with three components.
    GQCP::QCMatrix<double> h_x {2};
    // clang-format off
    h_x << 1.0, 2.0,
           3.0, 4.0;
    // clang-format on

    GQCP::QCMatrix<double> h_y {2};
    // clang-format off
    h_y << -1.0,  2.0,
            3.0, -4.0;
    // clang-format on

    GQCP::QCMatrix<double> h_z {2};
    // clang-format off
    h_z <<  0.0,  3.0,
            3.0,  0.0;
    // clang-format on

    GQCP::VectorRSQOneElectronOperator<double> h_op {{h_x, h_y, h_z}};


    // Define a vector to do the dot product with, and check the result.
    GQCP::VectorX<double> a {3};
    a << 1.0, 2.0, 3.0;

    GQCP::QCMatrix<double> dot_ref {2};
    // clang-format off
    dot_ref <<  -1.0,  15.0,
                18.0,  -4.0;
    // clang-format on

    BOOST_CHECK(h_op.dot(a).parameters().isApprox(dot_ref, 1.0e-12));
}
