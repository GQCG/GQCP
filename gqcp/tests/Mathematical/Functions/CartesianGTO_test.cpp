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

#define BOOST_TEST_MODULE "CartesianGTO"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Functions/CartesianGTO.hpp"


/**
 *  Check if the constructor throws as expected.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    const GQCP::CartesianExponents exponents {0, 0, 0};
    const GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();

    BOOST_CHECK_THROW(GQCP::CartesianGTO gto(-1.0, exponents, center), std::invalid_argument);  // exponent in the exponential cannot be negative
}


/**
 *  Check the evaluation of CartesianGTO with calculations by WolframAlpha.
 */
BOOST_AUTO_TEST_CASE(operator_call) {

    const GQCP::CartesianExponents exponents1 {1, 0, 1};
    const GQCP::Vector<double, 3> center1 {1.0, 0.0, -0.5};
    const GQCP::Vector<double, 3> r1 {0.0, 1.0, 0.0};
    const GQCP::CartesianGTO gto1 {1.0, exponents1, center1};

    BOOST_CHECK(std::abs(gto1(r1) - (-0.05269961228093216839)) < 1.0e-12);  // calculation by WolframAlpha


    const GQCP::CartesianExponents exponents2 {1, 2, 3};
    const GQCP::Vector<double, 3> center2 {0.0, 1.0, 0.0};
    const GQCP::Vector<double, 3> r2 {-1.0, -1.0, 1.0};
    const GQCP::CartesianGTO gto2 {2.5, exponents2, center2};

    BOOST_CHECK(std::abs(gto2(r2) - (-1.22360928200730315348e-06)) < 1.0e-12);  //calculation by WolframAlpha
}


/**
 *  Test the behaviour of CartesianGTO::operator==.
 */
BOOST_AUTO_TEST_CASE(operator_equals) {

    const double gaussian_exponent1 = 1.0;
    const double gaussian_exponent2 = 1.1;
    const GQCP::CartesianExponents cartesian_exponents1 {1, 0, 0};
    const GQCP::CartesianExponents cartesian_exponents2 {0, 1, 0};
    const GQCP::Vector<double, 3> center1 = GQCP::Vector<double, 3>::Zero(3);
    const GQCP::Vector<double, 3> center2 = GQCP::Vector<double, 3>::Ones(3);


    // Check for equality.
    BOOST_CHECK(GQCP::CartesianGTO(gaussian_exponent1, cartesian_exponents1, center1) == GQCP::CartesianGTO(gaussian_exponent1, cartesian_exponents1, center1));

    // Check for inequality if the Gaussian exponent doesn't match.
    BOOST_CHECK(!(GQCP::CartesianGTO(gaussian_exponent1, cartesian_exponents1, center1) == GQCP::CartesianGTO(gaussian_exponent2, cartesian_exponents1, center1)));

    // Check for inequality if the Cartesian exponents don't match.
    BOOST_CHECK(!(GQCP::CartesianGTO(gaussian_exponent1, cartesian_exponents1, center1) == GQCP::CartesianGTO(gaussian_exponent1, cartesian_exponents2, center1)));

    // Check for inequality if the center doesn't match.
    BOOST_CHECK(!(GQCP::CartesianGTO(gaussian_exponent1, cartesian_exponents1, center1) == GQCP::CartesianGTO(gaussian_exponent1, cartesian_exponents1, center2)));
}


/**
 *  Check the calculation of the CartesianGTO normalization factor with calculations by WolframAlpha.
 */
BOOST_AUTO_TEST_CASE(calculateNormalizationFactor) {

    const GQCP::CartesianExponents exponents1 {1, 0, 1};
    const GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();
    const GQCP::CartesianGTO gto1 {1.0, exponents1, center};

    BOOST_CHECK(std::abs(gto1.normalizationFactor() - 2.8508218814) < 1.0e-09);  // calculation by WolframAlpha

    const GQCP::CartesianExponents exponents2 {1, 2, 3};
    const GQCP::CartesianGTO gto2 {2.5, exponents2, center};

    BOOST_CHECK(std::abs(gto2.normalizationFactor() - 211.2315772257) < 1.0e-09);  // calculation by WolframAlpha
}


/**
 *  Check the calculation of the position derivative of a CartesianGTO.
 */
BOOST_AUTO_TEST_CASE(calculatePositionDerivative) {

    // First test case.
    const GQCP::CartesianExponents exponents1 {1, 0, 1};
    const GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();
    const GQCP::CartesianGTO gto1 {1.0, exponents1, center};

    // First test case - x-component
    const auto x_derivative1 = gto1.calculatePositionDerivative(GQCP::CartesianDirection::x);
    const double ref_coeff1_x1 = -2.0;
    const double ref_coeff2_x1 = 1.0;
    const GQCP::CartesianExponents ref_exp1_x1 {2, 0, 1};
    const GQCP::CartesianExponents ref_exp2_x1 {0, 0, 1};

    BOOST_CHECK(std::abs(x_derivative1.coefficients()[0] - ref_coeff1_x1) < 1.0e-12);
    BOOST_CHECK(std::abs(x_derivative1.coefficients()[1] - ref_coeff2_x1) < 1.0e-12);
    BOOST_CHECK(x_derivative1.functions()[0].cartesianExponents() == ref_exp1_x1);
    BOOST_CHECK(x_derivative1.functions()[1].cartesianExponents() == ref_exp2_x1);


    // First test case - y-component
    const auto y_derivative1 = gto1.calculatePositionDerivative(GQCP::CartesianDirection::y);
    const double ref_coeff1_y1 = -2.0;
    const GQCP::CartesianExponents ref_exp1_y1 {1, 1, 1};

    BOOST_CHECK(std::abs(y_derivative1.coefficients()[0] - ref_coeff1_y1) < 1.0e-12);
    BOOST_CHECK(y_derivative1.functions()[0].cartesianExponents() == ref_exp1_y1);


    // Second test case.
    const GQCP::CartesianExponents exponents2 {1, 2, 3};
    const GQCP::CartesianGTO gto2 {2.5, exponents2, center};

    // Second test case - x-component
    const auto x_derivative2 = gto2.calculatePositionDerivative(GQCP::CartesianDirection::x);
    const double ref_coeff1_x2 = -5.0;
    const double ref_coeff2_x2 = 1.0;
    const GQCP::CartesianExponents ref_exp1_x2 {2, 2, 3};
    const GQCP::CartesianExponents ref_exp2_x2 {0, 2, 3};

    BOOST_CHECK(std::abs(x_derivative2.coefficients()[0] - ref_coeff1_x2) < 1.0e-12);
    BOOST_CHECK(std::abs(x_derivative2.coefficients()[1] - ref_coeff2_x2) < 1.0e-12);
    BOOST_CHECK(x_derivative2.functions()[0].cartesianExponents() == ref_exp1_x2);
    BOOST_CHECK(x_derivative2.functions()[1].cartesianExponents() == ref_exp2_x2);

    // Second test case - z-component
    const auto z_derivative2 = gto2.calculatePositionDerivative(GQCP::CartesianDirection::z);
    const double ref_coeff1_z2 = -5.0;
    const double ref_coeff2_z2 = 3.0;
    const GQCP::CartesianExponents ref_exp1_z2 {1, 2, 4};
    const GQCP::CartesianExponents ref_exp2_z2 {1, 2, 2};

    BOOST_CHECK(std::abs(z_derivative2.coefficients()[0] - ref_coeff1_z2) < 1.0e-12);
    BOOST_CHECK(std::abs(z_derivative2.coefficients()[1] - ref_coeff2_z2) < 1.0e-12);
    BOOST_CHECK(z_derivative2.functions()[0].cartesianExponents() == ref_exp1_z2);
    BOOST_CHECK(z_derivative2.functions()[1].cartesianExponents() == ref_exp2_z2);
}
