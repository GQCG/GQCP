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
#define BOOST_TEST_MODULE "SQTwoElectronOperator"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "Utilities/linalg.hpp"
#include "Utilities/miscellaneous.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check the interface for constructing SQTwoElectronOperators from Tensors
 */
BOOST_AUTO_TEST_CASE(SQTwoElectronOperator_constructor) {

    // Check a correct constructor
    const GQCP::QCRankFourTensor<double> tensor {3};
    GQCP::ScalarSQTwoElectronOperator<double> O {tensor};


    // Check a faulty constructor
    GQCP::Tensor<double, 4> tensor2 {3, 3, 3, 2};
    BOOST_CHECK_THROW(GQCP::ScalarSQTwoElectronOperator<double> O2 {tensor2}, std::invalid_argument);
}


/**
 *  Check if the zero constructor really sets is parameters to all zeros
 */
BOOST_AUTO_TEST_CASE(SQTwoElectronOperator_zero_constructor) {

    const size_t dim = 2;
    GQCP::ScalarSQTwoElectronOperator<double> op {dim};

    // Create a reference zero tensor
    GQCP::QCRankFourTensor<double> ref {dim};

    BOOST_CHECK_EQUAL(op.dimension(), dim);
    BOOST_CHECK(op.parameters().isApprox(ref.setZero(), 1.0e-08));
}


/**
 *  Check if the formulas in effectiveOneElectronPartition are implemented correctly
 */
BOOST_AUTO_TEST_CASE(SQTwoElectronOperator_effectiveOneElectronPartition) {

    const size_t K = 4;
    auto K_ = static_cast<double>(K);

    // Set up toy 2-electron integrals
    GQCP::QCRankFourTensor<double> g_par {K};
    g_par.setZero();

    for (size_t i = 0; i < K; i++) {
        for (size_t j = 0; j < K; j++) {
            for (size_t k = 0; k < K; k++) {
                for (size_t l = 0; l < K; l++) {
                    g_par(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }

    GQCP::ScalarSQTwoElectronOperator<double> g {g_par};


    // Set up the reference effective one-electron integrals by manual calculation
    GQCP::QCMatrix<double> k_par_ref = GQCP::QCMatrix<double>::Zero(K, K);  // reference parameters
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            auto p_ = static_cast<double>(p) + 1;
            auto q_ = static_cast<double>(q) + 1;

            k_par_ref(p, q) = -K_ / 2 * (p_ + 8 * q_ + 3 * K_ + 3);
        }
    }


    BOOST_CHECK(k_par_ref.isApprox(g.effectiveOneElectronPartition().parameters(), 1.0e-08));
}


/**
 *  Check if calculateExpectationValue throws when necessary
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_throw) {

    const GQCP::ScalarSQTwoElectronOperator<double> g {2};

    const GQCP::TwoRDM<double> d_valid {2};
    const GQCP::TwoRDM<double> d_invalid {3};

    BOOST_CHECK_THROW(g.calculateExpectationValue(d_invalid), std::invalid_argument);
    BOOST_CHECK_NO_THROW(g.calculateExpectationValue(d_valid));
}


/**
 * Check whether or not calculateExpectationValue shows the correct behaviour
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_behaviour) {

    const size_t dim = 2;

    // Initialize a test tensor and convert it into an operator
    GQCP::QCRankFourTensor<double> T1 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }
    const GQCP::ScalarSQTwoElectronOperator<double> op(T1);

    // Initialize an alpha and beta density matrix, each one is chosen as a Hermitian matrix.
    GQCP::TwoRDM<double> d {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    d(i, j, k, l) = 1;
                }
            }
        }
    }

    // Initialize a reference value
    const double reference_expectation_value = 180.0;

    const auto expectation_value = op.calculateExpectationValue(d)(0);
    BOOST_CHECK(std::abs(expectation_value - reference_expectation_value) < 1.0e-08);
}


/**
 *  Check if addition of operators works as expected
 */
BOOST_AUTO_TEST_CASE(SQTwoElectronOperator_addition) {

    const size_t dim = 2;

    // Initialize two test tensors and convert them into operators
    GQCP::QCRankFourTensor<double> T1 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }
    const GQCP::ScalarSQTwoElectronOperator<double> op1 {T1};

    GQCP::QCRankFourTensor<double> T2 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T2(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }
    const GQCP::ScalarSQTwoElectronOperator<double> op2 {T2};


    // Initialize the reference and check the result
    GQCP::QCRankFourTensor<double> T_sum_ref {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T_sum_ref(i, j, k, l) = 2 * (i + 1) + 4 * (j + 1) + 8 * (k + 1) + 16 * (l + 1);
                }
            }
        }
    }

    const auto op_sum = op1 + op2;
    BOOST_CHECK(op_sum.parameters().isApprox(T_sum_ref, 1.0e-08));
}


/**
 *  Check if the scalar product with an operator works as expected
 */
BOOST_AUTO_TEST_CASE(SQTwoElectronOperator_scalar_product) {

    const size_t dim = 2;
    const double scalar = 2.0;

    // Initialize two test tensors and convert them into operators
    GQCP::QCRankFourTensor<double> T1 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }
    const GQCP::ScalarSQTwoElectronOperator<double> op1 {T1};

    const auto op_prod = scalar * op1;
    BOOST_CHECK(op_prod.parameters().isApprox((2 * T1), 1.0e-08));
}


/**
 *  Check if negating an operator works as expected
 */
BOOST_AUTO_TEST_CASE(SQTwoElectronOperator_negate) {

    const size_t dim = 2;

    // Initialize a test tensor and convert it into an operator
    GQCP::QCRankFourTensor<double> T1 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = 5;
                }
            }
        }
    }
    const GQCP::ScalarSQTwoElectronOperator<double> op1 {T1};

    // Initialize the reference and check the result
    GQCP::QCRankFourTensor<double> T_neg_ref {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T_neg_ref(i, j, k, l) = -5;
                }
            }
        }
    }

    const auto op_neg = -op1;
    BOOST_CHECK(op_neg.parameters().isApprox(T_neg_ref, 1.0e-08));
}


/**
 *  Check if taking the difference of two operators works as expected
 */
BOOST_AUTO_TEST_CASE(SQTwoElectronOperator_difference) {

    const size_t dim = 2;

    // Initialize two test tensors and convert them into operators
    GQCP::QCRankFourTensor<double> T1 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }
    const GQCP::ScalarSQTwoElectronOperator<double> op1 {T1};

    GQCP::QCRankFourTensor<double> T2 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T2(i, j, k, l) = 2 * (i + 1) + 4 * (j + 1) + 8 * (k + 1) + 16 * (l + 1);
                }
            }
        }
    }
    const GQCP::ScalarSQTwoElectronOperator<double> op2 {T2};

    const auto op_diff = op2 - op1;
    BOOST_CHECK(op_diff.parameters().isApprox(T1, 1.0e-08));
}


/**
 * Check whether or not the rotate with transformation matrix method works as expected
 */
BOOST_AUTO_TEST_CASE(rotate_with_unitary_transformation_matrix) {

    const size_t dim = 2;

    // Initialize a test tensor and convert it into an operator
    GQCP::QCRankFourTensor<double> T1 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }
    GQCP::ScalarSQTwoElectronOperator<double> op {T1};

    // Initialize a unitary transformation matrix
    GQCP::TransformationMatrix<double> U {dim};
    // clang-format off
    U << 1.0, 0.0,
         0.0, 1.0;
    // clang-format on

    op.rotate(U);
    BOOST_CHECK(op.parameters().isApprox(T1, 1.0e-08));
}


/**
 * Check whether or not the transform with transformation matrix method works as expected
 */
BOOST_AUTO_TEST_CASE(transform_with_transformation_matrix) {

    const size_t dim = 2;

    // Initialize a test tensor and convert it into an operator
    GQCP::QCRankFourTensor<double> T1 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = 1;
                }
            }
        }
    }
    GQCP::ScalarSQTwoElectronOperator<double> op {T1};

    // Initialize a transformation matrix
    GQCP::TransformationMatrix<double> T {dim};
    // clang-format off
    T << 2.0, 3.0,
         3.0, 4.0;
    // clang-format

    // Initialize a reference tensor
    GQCP::QCRankFourTensor<double> ref {dim};
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    if ((i + j + k + l) == 0) {
                        ref(i, j, k, l) = 625.0;
                    }
                    if ((i + j + k + l) == 1) {
                        ref(i, j, k, l) = 875.0;
                    }
                    if ((i + j + k + l) == 2) {
                        ref(i, j, k, l) = 1225.0;
                    }
                    if ((i + j + k + l) == 3) {
                        ref(i, j, k, l) = 1715.0;
                    }
                    if ((i + j + k + l) == 4) {
                        ref(i, j, k, l) = 2401.0;
                    }
                }
            }
        }
    }

    op.transform(T);
    BOOST_CHECK(op.parameters().isApprox(ref, 1.0e-08));
}


/**
 * Check whether or not the jacobi rotation method works as expected
 */
BOOST_AUTO_TEST_CASE(transform_with_jacobi_matrix) {

    const size_t dim = 2;

    // Initialize a test tensor and convert it into an operator
    GQCP::QCRankFourTensor<double> T1 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = 1;
                }
            }
        }
    }
    GQCP::ScalarSQTwoElectronOperator<double> op {T1};

    // Initialize a transformation matrix
    GQCP::JacobiRotationParameters J {1, 0, (boost::math::constants::pi<double>() / 2)};

    // Initialize a reference tensor
    GQCP::QCRankFourTensor<double> ref {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    if ((i + j + k + l) % 2 == 0) {
                        ref(i, j, k, l) = 1.0;
                    }
                    if ((i + j + k + l) % 2 != 0) {
                        ref(i, j, k, l) = -1.0;
                    }
                }
            }
        }
    }

    op.rotate(J);
    BOOST_CHECK(op.parameters().isApprox(ref, 1.0e-08));
}
