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
#define BOOST_TEST_MODULE "OneElectronOperator"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


#include "Operator/OneElectronOperator.hpp"

#include "CartesianGTO.hpp"
#include "JacobiRotationParameters.hpp"
#include "utilities/miscellaneous.hpp"


BOOST_AUTO_TEST_CASE ( OneElectronOperator_constructor ) {

    // Check a correct constructor
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(4, 4);
    GQCP::OneElectronOperator<double> O (matrix);


    // Check a faulty constructor
    Eigen::MatrixXd matrix2 = Eigen::MatrixXd::Zero(3, 4);
    BOOST_CHECK_THROW(GQCP::OneElectronOperator<double> O2 (matrix2), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( operator_plus ) {
    
    // Construct two OneElectronOperators
    size_t K = 5;
    Eigen::MatrixXd matrix1 = Eigen::MatrixXd::Random(K, K);
    Eigen::MatrixXd matrix2 = Eigen::MatrixXd::Random(K, K);
    
    GQCP::OneElectronOperator<double> M1 (matrix1);
    GQCP::OneElectronOperator<double> M2 (matrix2);
    
    
    // Check if operator+ works
    BOOST_CHECK((M1 + M2).isApprox(matrix1 + matrix2, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_transform_trivial ) {

    // Let's test a trivial transformation: i.e. with T being a unit matrix
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    GQCP::OneElectronOperator<double> H (h);

    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    H.transform(GQCP::SquareMatrix<double>(T));

    BOOST_CHECK(H.isApprox(h, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_transform_and_inverse ) {

    // Let's test if, if we transform h with T and then with T_inverse, we get effectively do nothing
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    GQCP::OneElectronOperator<double> H (h);

    Eigen::MatrixXd T (3, 3);
    T << 1,  0,  0,
         0, -2,  0,
         0,  0,  3;
    Eigen::MatrixXd T_inverse = T.inverse();


    H.transform(GQCP::SquareMatrix<double>(T));
    H.transform(GQCP::SquareMatrix<double>(T_inverse));

    BOOST_CHECK(H.isApprox(h, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_rotate_throws ) {

    // Create a random OneElectronOperator
    size_t dim = 3;
    GQCP::OneElectronOperator<double> M (Eigen::MatrixXd::Random(dim, dim));


    // Check if a non-unitary matrix as transformation matrix causes a throw
    Eigen::MatrixXd U (Eigen::MatrixXd::Random(dim, dim));
    BOOST_CHECK_THROW(M.rotate(GQCP::SquareMatrix<double>(U)), std::invalid_argument);


    // Check if a unitary matrix as transformation matrix is accepted
    M.rotate(GQCP::SquareMatrix<double>(GQCP::Matrix<double>::Identity(dim, dim)));
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_rotate_JacobiRotationParameters ) {

    // Create a random OneElectronOperator
    size_t dim = 5;
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(dim, dim);
    GQCP::OneElectronOperator<double> M1 (m);
    GQCP::OneElectronOperator<double> M2 (m);


    // Check that using a Jacobi transformation (rotation) matrix as U is equal to the custom transformation (rotation)
    // with custom JacobiRotationParameters
    GQCP::JacobiRotationParameters jacobi_rotation_parameters (4, 2, 56.81);

    auto U = GQCP::jacobiRotationMatrix(jacobi_rotation_parameters, dim);


    M1.rotate(jacobi_rotation_parameters);
    M2.rotate(GQCP::SquareMatrix<double>(U));


    BOOST_CHECK(M1.isApprox(M2, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( transformation_GTO_integrals ) {

    // Test the transformation of a one-electron operator consisting of GTOs

    // Build up the one-electron operator with linear combinations of GTOs
    Eigen::Vector3d center = Eigen::Vector3d::Zero();

    double coeff1 = 1.0;
    std::array<size_t, 3> exp1 {1, 0, 0};
    GQCP::CartesianGTO gto1 (1.0, exp1, center);

    double coeff2 = 2.0;
    std::array<size_t, 3> exp2 {0, 1, 0};
    GQCP::CartesianGTO gto2 (2.0, exp2, center);

    double coeff3 = -1.0;
    std::array<size_t, 3> exp3 {1, 1, 0};
    GQCP::CartesianGTO gto3 (1.0, exp3, center);

    double coeff4 = 1.0;
    std::array<size_t, 3> exp4 {0, 0, 2};
    GQCP::CartesianGTO gto4 (3.0, exp4, center);

    double coeff5 = 2.5;
    std::array<size_t, 3> exp5 {1, 1, 1};
    GQCP::CartesianGTO gto5 (0.5, exp5, center);

    double coeff6 = -1.0;
    std::array<size_t, 3> exp6 {0, 1, 1};
    GQCP::CartesianGTO gto6 (0.5, exp6, center);


    std::vector<double> lc1_coeff {coeff1};
    std::vector<double> lc2_coeff {coeff2, coeff3};
    std::vector<double> lc3_coeff {coeff4, coeff5};
    std::vector<double> lc4_coeff {coeff6};

    std::vector<GQCP::CartesianGTO> lc1_funcs {gto1};
    std::vector<GQCP::CartesianGTO> lc2_funcs {gto2, gto3};
    std::vector<GQCP::CartesianGTO> lc3_funcs {gto4, gto5};
    std::vector<GQCP::CartesianGTO> lc4_funcs {gto6};


    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc1 (lc1_coeff, lc1_funcs);
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc2 (lc2_coeff, lc2_funcs);
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc3 (lc3_coeff, lc3_funcs);
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc4 (lc4_coeff, lc4_funcs);


    Eigen::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho;
    rho << lc1, lc2, lc3, lc4;


    // Create the transformation matrix
    Eigen::Matrix2d T = Eigen::Matrix2d::Zero();
    T << 2.0, 1.0, 1.0, 0.0;


    Eigen::Matrix<GQCP::LinearCombination<double, GQCP::CartesianGTO>, 2, 2> rho_transformed = T.adjoint() * rho * T;


    // Check the coefficients of the results
    std::vector<double> ref_coeff_result_01 {2.0, 1.0, 2.5};
    auto coeff_result_01 = rho_transformed(0,1).get_coefficients();
    for (size_t i = 0; i < ref_coeff_result_01.size(); i++) {
        BOOST_CHECK(std::abs(ref_coeff_result_01[i] - coeff_result_01[i]) < 1.0e-12);
    }
    BOOST_CHECK(ref_coeff_result_01.size() == coeff_result_01.size());

    std::vector<double> ref_coeff_result_11 {1.0};
    auto coeff_result_11 = rho_transformed(1,1).get_coefficients();
    for (size_t i = 0; i < ref_coeff_result_11.size(); i++) {
        BOOST_CHECK(std::abs(ref_coeff_result_11[i] - coeff_result_11[i]) < 1.0e-12);
    }
    BOOST_CHECK(ref_coeff_result_11.size() == coeff_result_11.size());
}
