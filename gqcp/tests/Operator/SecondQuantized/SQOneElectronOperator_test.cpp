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

#define BOOST_TEST_MODULE "SQOneElectronOperator"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Utilities/miscellaneous.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check if the transformation of a one-electron operator consisting of linear combinations of GTOs can be supported through the underlying scalar types
 */
BOOST_AUTO_TEST_CASE(SQOneElectronOperator_of_GTOs_transform) {

    // Create a toy operator of linear combinations of GTOs that correspond to a manual calculation.
    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();

    const double coeff1 = 1.0;
    const GQCP::CartesianGTO gto1 {1.0, {1, 0, 0}, center};

    const double coeff2 = 2.0;
    const GQCP::CartesianGTO gto2 {2.0, {0, 1, 0}, center};

    const double coeff3 = -1.0;
    const GQCP::CartesianGTO gto3 {1.0, {1, 1, 0}, center};

    const double coeff4 = 1.0;
    const GQCP::CartesianGTO gto4 {3.0, {0, 0, 2}, center};

    const double coeff5 = 2.5;
    const GQCP::CartesianGTO gto5 {0.5, {1, 1, 1}, center};

    const double coeff6 = -1.0;
    const GQCP::CartesianGTO gto6 {2.5, {0, 1, 1}, center};


    const GQCP::LinearCombination<double, GQCP::CartesianGTO> lc1 {coeff1, gto1};
    const GQCP::LinearCombination<double, GQCP::CartesianGTO> lc2 {{coeff2, coeff3}, {gto2, gto3}};
    const GQCP::LinearCombination<double, GQCP::CartesianGTO> lc3 {{coeff4, coeff5}, {gto4, gto5}};
    const GQCP::LinearCombination<double, GQCP::CartesianGTO> lc4 {coeff6, gto6};


    GQCP::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho_par;
    // clang-format off
    rho_par << lc1, lc2,
               lc3, lc4;
    // clang-format on


    // Create the transformation matrix and transform the operator manually: .transform() for does not work yet
    GQCP::Matrix<double, 2, 2> T = GQCP::Matrix<double, 2, 2>::Zero();
    // clang-format off
    T << 2.0, 1.0,
         1.0, 0.0;
    // clang-format on

    const GQCP::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho_transformed_par = T.adjoint() * rho_par * T;


    // Check the coefficients of the transformed operator with a manual calculation
    const std::vector<double> ref_coeff_result_01 {2.0, 1.0, 2.5};
    const auto coeff_result_01 = rho_transformed_par(0, 1).coefficients();
    for (size_t i = 0; i < ref_coeff_result_01.size(); i++) {
        BOOST_CHECK(std::abs(ref_coeff_result_01[i] - coeff_result_01[i]) < 1.0e-12);
    }

    BOOST_CHECK(ref_coeff_result_01.size() == coeff_result_01.size());

    const std::vector<double> ref_coeff_result_11 {1.0};
    const auto coeff_result_11 = rho_transformed_par(1, 1).coefficients();
    for (size_t i = 0; i < ref_coeff_result_11.size(); i++) {
        BOOST_CHECK(std::abs(ref_coeff_result_11[i] - coeff_result_11[i]) < 1.0e-12);
    }
    BOOST_CHECK(ref_coeff_result_11.size() == coeff_result_11.size());
}


/**
 *  Check if we can evaluate a SQOneElectronOperator consisting of GTOs in a given point r
 */

BOOST_AUTO_TEST_CASE(SQOneElectronOperator_of_GTOs_evaluate) {

    // Create a toy operator of linear combinations of GTOs that correspond to a manual calculation.
    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();

    const double coeff1 = 1.0;
    const GQCP::CartesianGTO gto1 {1.0, {1, 0, 0}, center};

    const double coeff2 = 2.0;
    const GQCP::CartesianGTO gto2 {2.0, {0, 1, 0}, center};

    const double coeff3 = -1.0;
    const GQCP::CartesianGTO gto3 {1.0, {1, 1, 0}, center};

    const double coeff4 = 1.0;
    const GQCP::CartesianGTO gto4 {3.0, {0, 0, 2}, center};

    const double coeff5 = 2.5;
    const GQCP::CartesianGTO gto5 {0.5, {1, 1, 1}, center};

    const double coeff6 = -1.0;
    const GQCP::CartesianGTO gto6 {2.5, {0, 1, 1}, center};


    const GQCP::LinearCombination<double, GQCP::CartesianGTO> lc1 {coeff1, gto1};
    const GQCP::LinearCombination<double, GQCP::CartesianGTO> lc2 {{coeff2, coeff3}, {gto2, gto3}};
    const GQCP::LinearCombination<double, GQCP::CartesianGTO> lc3 {{coeff4, coeff5}, {gto4, gto5}};
    const GQCP::LinearCombination<double, GQCP::CartesianGTO> lc4 {coeff6, gto6};


    GQCP::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho_par;
    // clang-format off
    rho_par << lc1, lc2,
               lc3, lc4;
    // clang-format on


    // Create the transformation matrix and transform the operator manually: .transform() for does not work yet
    GQCP::Matrix<double, 2, 2> T = GQCP::Matrix<double, 2, 2>::Zero();
    // clang-format off
    T << 2.0, 1.0,
         1.0, 0.0;
    // clang-format on

    const GQCP::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho_transformed_par = T.adjoint() * rho_par * T;
    const GQCP::ScalarSQOneElectronOperator<GQCP::LinearCombination<double, GQCP::CartesianGTO>> rho {rho_transformed_par};


    // Evaluate the operator of GTOs at the given point r
    const GQCP::Vector<double, 3> r {1.0, 1.0, 1.0};
    const auto rho_evaluated_par = rho.evaluate(r).parameters();


    // Read in the reference solution and check the results.
    GQCP::SquareMatrix<double> ref_rho_evaluated_par = GQCP::SquareMatrix<double>::Zero(2);
    double ref_rho_evaluated_00 = 4 * std::exp(-3.0) + 4 * std::exp(-6.0) - 2 * std::exp(-3.0) + 2 * std::exp(-9.0) + 5 * std::exp(-1.5) - 1 * std::exp(-7.5);
    double ref_rho_evaluated_01 = 2 * std::exp(-3.0) + 1 * std::exp(-9.0) + 2.5 * std::exp(-1.5);
    double ref_rho_evaluated_10 = 2 * std::exp(-3.0) + 2 * std::exp(-6.0) - 1 * std::exp(-3.0);
    double ref_rho_evaluated_11 = 1 * std::exp(-3.0);
    // clang-format off
    ref_rho_evaluated_par << ref_rho_evaluated_00, ref_rho_evaluated_01, 
                             ref_rho_evaluated_10, ref_rho_evaluated_11;
    // clang-format on

    BOOST_CHECK(ref_rho_evaluated_par.isApprox(rho_evaluated_par, 1.0e-12));
}
