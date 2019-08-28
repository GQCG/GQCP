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
#define BOOST_TEST_MODULE "SQOneElectronOperator"

#include <boost/test/unit_test.hpp>

#include "Basis/CartesianGTO.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Utilities/miscellaneous.hpp"


/**
 *  Check the construction of one-electron operators from matrices
 */
BOOST_AUTO_TEST_CASE ( SQOneElectronOperator_constructor ) {

    // Check a correct constructor
    const auto square_matrix = GQCP::SquareMatrix<double>::Zero(4, 4);
    GQCP::ScalarSQOneElectronOperator<double> O ({square_matrix});


    // Check a faulty constructor
    GQCP::MatrixX<double> matrix = GQCP::MatrixX<double>::Zero(3, 4);
    BOOST_CHECK_THROW(GQCP::ScalarSQOneElectronOperator<double> O2 ({matrix}), std::invalid_argument);
}


/**
 *  Check if the transformation of a one-electron operator consisting of linear combinations of GTOs can be supported through the underlying scalar types
 */
BOOST_AUTO_TEST_CASE ( SQOneElectronOperator_of_GTOs_transform ) {

    // Create a toy operator of linear combinations of GTOs that correspond to a manual calculation
    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();

    double coeff1 = 1.0;
    GQCP::CartesianGTO gto1 (1.0, {1, 0, 0}, center);
    double N1 = gto1.calculateNormalizationFactor();

    double coeff2 = 2.0;
    GQCP::CartesianGTO gto2 (2.0, {0, 1, 0}, center);
    double N2 = gto2.calculateNormalizationFactor();

    double coeff3 = -1.0;
    GQCP::CartesianGTO gto3 (1.0, {1, 1, 0}, center);
    double N3 = gto3.calculateNormalizationFactor();

    double coeff4 = 1.0;
    GQCP::CartesianGTO gto4 (3.0, {0, 0, 2}, center);
    double N4 = gto4.calculateNormalizationFactor();

    double coeff5 = 2.5;
    GQCP::CartesianGTO gto5 (0.5, {1, 1, 1}, center);
    double N5 = gto5.calculateNormalizationFactor();

    double coeff6 = -1.0;
    GQCP::CartesianGTO gto6 (2.5, {0, 1, 1}, center);
    double N6 = gto6.calculateNormalizationFactor();


    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc1 (coeff1, gto1);
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc2 ({coeff2, coeff3}, {gto2, gto3});
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc3 ({coeff4, coeff5}, {gto4, gto5});
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc4 (coeff6, gto6);


    GQCP::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho_par;
    rho_par << lc1, lc2, 
               lc3, lc4;


    // Create the transformation matrix and transform the operator manually: .transform() for does not work yet
    GQCP::Matrix<double, 2, 2> T = GQCP::Matrix<double, 2, 2>::Zero();
    T << 2.0, 1.0, 
         1.0, 0.0;

    Eigen::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho_transformed_par = T.adjoint() * rho_par * T;


    // Check the coefficients of the transformed operator with a manual calculation
    const std::vector<double> ref_coeff_result_01 {2.0, 1.0, 2.5};
    auto coeff_result_01 = rho_transformed_par(0,1).get_coefficients();
    for (size_t i = 0; i < ref_coeff_result_01.size(); i++) {
        BOOST_CHECK(std::abs(ref_coeff_result_01[i] - coeff_result_01[i]) < 1.0e-12);
    }
    BOOST_CHECK(ref_coeff_result_01.size() == coeff_result_01.size());

    const std::vector<double> ref_coeff_result_11 {1.0};
    auto coeff_result_11 = rho_transformed_par(1,1).get_coefficients();
    for (size_t i = 0; i < ref_coeff_result_11.size(); i++) {
        BOOST_CHECK(std::abs(ref_coeff_result_11[i] - coeff_result_11[i]) < 1.0e-12);
    }
    BOOST_CHECK(ref_coeff_result_11.size() == coeff_result_11.size());
}


/**
 *  Check if we can evaluate a SQOneElectronOperator consisting of GTOs in a given point r
 */

BOOST_AUTO_TEST_CASE ( SQOneElectronOperator_of_GTOs_evaluate ) {

    // Create a toy operator of linear combinations of GTOs that correspond to a manual calculation. This is a copy-paste of the previous test: I can't make an easy auxiliary function because I also need N1, N2, ...
    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();

    double coeff1 = 1.0;
    GQCP::CartesianGTO gto1 (1.0, {1, 0, 0}, center);
    double N1 = gto1.calculateNormalizationFactor();

    double coeff2 = 2.0;
    GQCP::CartesianGTO gto2 (2.0, {0, 1, 0}, center);
    double N2 = gto2.calculateNormalizationFactor();

    double coeff3 = -1.0;
    GQCP::CartesianGTO gto3 (1.0, {1, 1, 0}, center);
    double N3 = gto3.calculateNormalizationFactor();

    double coeff4 = 1.0;
    GQCP::CartesianGTO gto4 (3.0, {0, 0, 2}, center);
    double N4 = gto4.calculateNormalizationFactor();

    double coeff5 = 2.5;
    GQCP::CartesianGTO gto5 (0.5, {1, 1, 1}, center);
    double N5 = gto5.calculateNormalizationFactor();

    double coeff6 = -1.0;
    GQCP::CartesianGTO gto6 (2.5, {0, 1, 1}, center);
    double N6 = gto6.calculateNormalizationFactor();


    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc1 (coeff1, gto1);
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc2 ({coeff2, coeff3}, {gto2, gto3});
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc3 ({coeff4, coeff5}, {gto4, gto5});
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc4 (coeff6, gto6);


    GQCP::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho_par;
    rho_par << lc1, lc2, 
               lc3, lc4;


    // Create the transformation matrix and transform the operator manually: .transform() for does not work yet
    GQCP::Matrix<double, 2, 2> T = GQCP::Matrix<double, 2, 2>::Zero();
    T << 2.0, 1.0, 
         1.0, 0.0;

    Eigen::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho_transformed_par = T.adjoint() * rho_par * T;
    GQCP::ScalarSQOneElectronOperator<GQCP::LinearCombination<double, GQCP::CartesianGTO>> rho ({rho_transformed_par});


    // Evaluate the operator of GTOs at the given point r
    GQCP::Vector<double, 3> r = GQCP::Vector<double, 3>::Zero();
    r << 1.0, 1.0, 1.0;

    const auto rho_evaluated_par = rho.evaluate(r).parameters();


    // Read in the reference solution and check the results
    GQCP::ChemicalMatrix<double> ref_rho_evaluated_par = GQCP::ChemicalMatrix<double>::Zero(2, 2);
    double ref_rho_evaluated_00 = 4*N1 * std::exp(-3.0) + 4*N2 * std::exp(-6.0) - 2*N3* std::exp(-3.0) + 2*N4 * std::exp(-9.0) + 5*N5 * std::exp(-1.5) - 1*N6 * std::exp(-7.5);
    double ref_rho_evaluated_01 = 2*N1 * std::exp(-3.0) + 1*N4 * std::exp(-9.0) + 2.5*N5 * std::exp(-1.5);
    double ref_rho_evaluated_10 = 2*N1 * std::exp(-3.0) + 2*N2 * std::exp(-6.0) - 1*N3 * std::exp(-3.0);
    double ref_rho_evaluated_11 = 1*N1 * std::exp(-3.0);
    ref_rho_evaluated_par << ref_rho_evaluated_00, ref_rho_evaluated_01, 
                             ref_rho_evaluated_10, ref_rho_evaluated_11;

    BOOST_CHECK(ref_rho_evaluated_par.isApprox(rho_evaluated_par, 1.0e-12));
}
