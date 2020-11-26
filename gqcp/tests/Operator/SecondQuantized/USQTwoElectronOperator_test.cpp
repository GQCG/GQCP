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

#define BOOST_TEST_MODULE "USQTwoElectronOperator"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/USQTwoElectronOperator.hpp"
#include "Utilities/miscellaneous.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check the interface for constructing USQTwoElectronOperators from Tensors
 */
BOOST_AUTO_TEST_CASE(USQTwoElectronOperator_constructor) {

    // Check a correct constructor.
    const GQCP::SquareRankFourTensor<double> tensor {3};
    GQCP::ScalarUSQTwoElectronOperator<double> O {tensor, tensor, tensor, tensor};  // All dimensions are equal.


    // Check a faulty constructor.
    GQCP::Tensor<double, 4> tensor2 {3, 3, 3, 2};
    BOOST_CHECK_THROW(GQCP::ScalarUSQTwoElectronOperator<double> O2(tensor2, tensor2, tensor2, tensor2), std::invalid_argument);  // All tensors have the wrong dimension.
    BOOST_CHECK_THROW(GQCP::ScalarUSQTwoElectronOperator<double> O2(tensor, tensor2, tensor, tensor), std::invalid_argument);     // Tensor 2 has the wrong dimension.
    BOOST_CHECK_THROW(GQCP::ScalarUSQTwoElectronOperator<double> O2(tensor2, tensor, tensor2, tensor), std::invalid_argument);    // Tensor 1 and 3 have the wrong dimension.
    BOOST_CHECK_THROW(GQCP::ScalarUSQTwoElectronOperator<double> O2(tensor2, tensor2, tensor2, tensor), std::invalid_argument);   // Tensor 1, 2, and 3 have the wrong dimension.
}


/**
 *  Check if the zero named constructor really sets its parameters to all zeros.
 */
BOOST_AUTO_TEST_CASE(USQTwoElectronOperator_zero_constructor) {

    const size_t K = 2;
    const auto op = GQCP::ScalarUSQTwoElectronOperator<double>::Zero(K);  // Should initialize to zeros.

    // Create a reference zero tensor.
    GQCP::SquareRankFourTensor<double> zero {K};
    zero.setZero();

    // Check if the number of orbitals is correct.
    BOOST_CHECK_EQUAL(op.numberOfOrbitals(GQCP::Spin::alpha, GQCP::Spin::alpha), K);
    BOOST_CHECK_EQUAL(op.numberOfOrbitals(GQCP::Spin::alpha, GQCP::Spin::beta), K);
    BOOST_CHECK_EQUAL(op.numberOfOrbitals(GQCP::Spin::beta, GQCP::Spin::alpha), K);
    BOOST_CHECK_EQUAL(op.numberOfOrbitals(GQCP::Spin::beta, GQCP::Spin::beta), K);


    BOOST_CHECK(op.alphaAlpha().parameters().isApprox(zero, 1.0e-08));
    BOOST_CHECK(op.alphaBeta().parameters().isApprox(zero, 1.0e-08));
    BOOST_CHECK(op.betaAlpha().parameters().isApprox(zero, 1.0e-08));
    BOOST_CHECK(op.betaBeta().parameters().isApprox(zero, 1.0e-08));
}


/**
 *  Check if calculateExpectationValue throws when necessary.
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_throw) {

    const size_t dim = 2;

    // Initialize a test operator.
    const auto g = GQCP::ScalarUSQTwoElectronOperator<double>::Zero(dim);

    // Initialize valid and invalid test 2-DMs.
    const GQCP::PureSpinResolved2DMComponent<double> d_valid_pure {2};
    const GQCP::PureSpinResolved2DMComponent<double> d_invalid_pure {3};
    const GQCP::MixedSpinResolved2DMComponent<double> d_valid_mixed {2};
    const GQCP::MixedSpinResolved2DMComponent<double> d_invalid_mixed {3};

    const GQCP::SpinResolved2DM<double> d_invalid_aa {d_invalid_pure, d_valid_mixed, d_valid_mixed, d_valid_pure};
    const GQCP::SpinResolved2DM<double> d_invalid_ab {d_valid_pure, d_invalid_mixed, d_valid_mixed, d_valid_pure};
    const GQCP::SpinResolved2DM<double> d_invalid_ba {d_valid_pure, d_valid_mixed, d_invalid_mixed, d_valid_pure};
    const GQCP::SpinResolved2DM<double> d_invalid_bb {d_valid_pure, d_valid_mixed, d_valid_mixed, d_invalid_pure};

    const GQCP::SpinResolved2DM<double> d_valid_spin_resolved {d_valid_pure, d_valid_mixed, d_valid_mixed, d_valid_pure};

    // Check if the calculateExpectationValue calls throw when expected.
    BOOST_CHECK_THROW(g.calculateExpectationValue(d_invalid_aa), std::invalid_argument);  // Wrong fist component.
    BOOST_CHECK_THROW(g.calculateExpectationValue(d_invalid_ab), std::invalid_argument);  // Wrong second component.
    BOOST_CHECK_THROW(g.calculateExpectationValue(d_invalid_ba), std::invalid_argument);  // Wrong third component.
    BOOST_CHECK_THROW(g.calculateExpectationValue(d_invalid_bb), std::invalid_argument);  // Wrong fourth component.

    BOOST_CHECK_NO_THROW(g.calculateExpectationValue(d_valid_spin_resolved));
}


/**
 * Check whether or not calculateExpectationValue shows the correct behaviour.
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_behaviour) {

    const size_t dim = 2;

    // Initialize test two-electron operators.
    auto T_pure = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T_pure(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }

    auto T_mixed = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T_mixed(i, j, k, l) = (i + 1) + (j + 1) + (k + 1) + (l + 1);
                }
            }
        }
    }
    const GQCP::ScalarUSQTwoElectronOperator<double> op {T_pure, T_mixed, T_mixed, T_pure};


    // Initialize test 2-DMs: each one is chosen to have the correct four-index symmetries.
    auto d1 = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    d1(i, j, k, l) = 1;
                }
            }
        }
    }

    auto d2 = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    d2(i, j, k, l) = 2;
                }
            }
        }
    }
    const GQCP::SpinResolved2DM<double> d {d1, d1, d2, d2};


    // Initialize a reference value and check the result.
    const double reference_expectation_value = 684.0;

    const double expectation_value = op.calculateExpectationValue(d);  // A scalar-StorageArray can be implicitly converted to its underlying scalar.
    BOOST_CHECK(std::abs(expectation_value - reference_expectation_value) < 1.0e-08);
}


/**
 *  Check if addition of unrestricted two-electron operators works as expected.
 */
BOOST_AUTO_TEST_CASE(USQTwoElectronOperator_addition) {

    const size_t dim = 2;

    // Initialize a test unrestricted two-electron operator.
    auto T1 = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }
    const GQCP::ScalarUSQTwoElectronOperator<double> op1 {T1, T1, T1, T1};


    // Initialize the reference and check the result.
    auto T_sum_ref = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T_sum_ref(i, j, k, l) = 2 * (i + 1) + 4 * (j + 1) + 8 * (k + 1) + 16 * (l + 1);
                }
            }
        }
    }

    auto op_sum = op1 + op1;
    BOOST_CHECK(op_sum.alphaAlpha().parameters().isApprox(T_sum_ref, 1.0e-08));
    BOOST_CHECK(op_sum.alphaBeta().parameters().isApprox(T_sum_ref, 1.0e-08));
    BOOST_CHECK(op_sum.betaAlpha().parameters().isApprox(T_sum_ref, 1.0e-08));
    BOOST_CHECK(op_sum.betaBeta().parameters().isApprox(T_sum_ref, 1.0e-08));
}


/**
 *  Check if the scalar multiplication with an unrestricted two-electron operator works as expected.
 */
BOOST_AUTO_TEST_CASE(USQTwoElectronOperator_scalar_product) {

    const size_t dim = 2;

    // Initialize a test unrestricted two-electron operator.
    GQCP::SquareRankFourTensor<double> T1 {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }
    const GQCP::ScalarUSQTwoElectronOperator<double> op1 {T1, T1, T1, T1};


    // Initialize the reference and check the result.
    auto T_mult_ref = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T_mult_ref(i, j, k, l) = 2 * (i + 1) + 4 * (j + 1) + 8 * (k + 1) + 16 * (l + 1);
                }
            }
        }
    }

    auto op_prod = 2.0 * op1;
    BOOST_CHECK(op_prod.alphaAlpha().parameters().isApprox(T_mult_ref, 1.0e-08));
    BOOST_CHECK(op_prod.alphaBeta().parameters().isApprox(T_mult_ref, 1.0e-08));
    BOOST_CHECK(op_prod.betaAlpha().parameters().isApprox(T_mult_ref, 1.0e-08));
    BOOST_CHECK(op_prod.betaBeta().parameters().isApprox(T_mult_ref, 1.0e-08));
}


/**
 *  Check whether or not rotating with a trivial unitary transformation has no apparent effect.
 */
BOOST_AUTO_TEST_CASE(rotate_with_unitary_transformation) {

    const size_t dim = 2;

    // Initialize a test tensor and convert it into an operator.
    auto T1 = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }
    GQCP::ScalarUSQTwoElectronOperator<double> op {T1, T1, T1, T1};

    // Initialize the identity unitary transformation and use it for the rotation.
    const auto U = GQCP::UTransformation<double>::Identity(dim);

    op.rotate(U);
    BOOST_CHECK(op.alphaAlpha().parameters().isApprox(T1, 1.0e-08));
    BOOST_CHECK(op.alphaBeta().parameters().isApprox(T1, 1.0e-08));
    BOOST_CHECK(op.betaAlpha().parameters().isApprox(T1, 1.0e-08));
    BOOST_CHECK(op.betaBeta().parameters().isApprox(T1, 1.0e-08));
}


/**
 * Check whether or not the transformation method works as expected.
 */
BOOST_AUTO_TEST_CASE(transform) {

    const size_t dim = 2;

    // Initialize a test tensor and convert it into an operator.
    auto T1 = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = 1;
                }
            }
        }
    }
    GQCP::ScalarUSQTwoElectronOperator<double> op {T1, T1, T1, T1};

    // Initialize a test transformation, equal for alpha and beta.
    GQCP::SquareMatrix<double> T_component {dim};
    // clang-format off
    T_component << 2.0, 3.0,
                   3.0, 4.0;
    // clang-format on
    const auto T = GQCP::UTransformation<double>::FromEqual(T_component);


    // Initialize a reference tensor, transform the original operator and check the result.
    auto ref = GQCP::SquareRankFourTensor<double>::Zero(dim);
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
    BOOST_CHECK(op.alphaAlpha().parameters().isApprox(ref, 1.0e-08));
    BOOST_CHECK(op.alphaBeta().parameters().isApprox(ref, 1.0e-08));
    BOOST_CHECK(op.betaAlpha().parameters().isApprox(ref, 1.0e-08));
    BOOST_CHECK(op.betaBeta().parameters().isApprox(ref, 1.0e-08));
}


/**
 * Check whether or not the jacobi rotation method works as expected.
 */
BOOST_AUTO_TEST_CASE(transform_with_jacobi_matrix) {

    const size_t dim = 2;

    // Initialize a test tensor and convert it into an operator.
    auto T1 = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    T1(i, j, k, l) = 1;
                }
            }
        }
    }
    GQCP::ScalarUSQTwoElectronOperator<double> op {T1, T1, T1, T1};

    // Initialize a Jacobi rotation.
    const GQCP::JacobiRotation J_component {1, 0, (boost::math::constants::pi<double>() / 2)};
    const auto J = GQCP::UJacobiRotation::FromEqual(J_component);

    // Initialize a reference tensor.
    auto ref = GQCP::SquareRankFourTensor<double>::Zero(dim);

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
    BOOST_CHECK(op.alphaAlpha().parameters().isApprox(ref, 1.0e-08));
    BOOST_CHECK(op.alphaBeta().parameters().isApprox(ref, 1.0e-08));
    BOOST_CHECK(op.betaAlpha().parameters().isApprox(ref, 1.0e-08));
    BOOST_CHECK(op.betaBeta().parameters().isApprox(ref, 1.0e-08));
}
