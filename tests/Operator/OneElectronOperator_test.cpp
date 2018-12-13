// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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


#include "Operator/OneElectronOperator.hpp"

#include "JacobiRotationParameters.hpp"
#include "utilities/miscellaneous.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( OneElectronOperator_constructor ) {

    // Check a correct constructor
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(4, 4);
    GQCP::OneElectronOperator O (matrix);


    // Check a faulty constructor
    Eigen::MatrixXd matrix2 = Eigen::MatrixXd::Zero(3, 4);
    BOOST_CHECK_THROW(GQCP::OneElectronOperator O2 (matrix2), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( operator_plus ) {
    
    // Construct two OneElectronOperators
    size_t K = 5;
    Eigen::MatrixXd matrix1 = Eigen::MatrixXd::Random(K, K);
    Eigen::MatrixXd matrix2 = Eigen::MatrixXd::Random(K, K);
    
    GQCP::OneElectronOperator M1 (matrix1);
    GQCP::OneElectronOperator M2 (matrix2);
    
    
    // Check if operator+ works
    BOOST_CHECK((M1 + M2).get_matrix_representation().isApprox(matrix1 + matrix2, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( isEqualTo ) {

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3, 3);
    GQCP::OneElectronOperator O1 (A);
    GQCP::OneElectronOperator O2 (A);
    BOOST_CHECK(O1.isEqualTo(O2, 1.0e-05));

    GQCP::OneElectronOperator O3 (2*A);
    BOOST_CHECK(!(O3.isEqualTo(O1)));
}


BOOST_AUTO_TEST_CASE ( operator_equals ) {

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3, 3);
    GQCP::OneElectronOperator O1 (A);
    GQCP::OneElectronOperator O2 (A);
    BOOST_CHECK(O1 == O2);
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_transform_trivial ) {

    // Let's test a trivial transformation: i.e. with T being a unit matrix
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    GQCP::OneElectronOperator H (h);

    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    H.transform(T);

    BOOST_CHECK(H.get_matrix_representation().isApprox(h, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_transform_and_inverse ) {

    // Let's test if, if we transform h with T and then with T_inverse, we get effectively do nothing
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    GQCP::OneElectronOperator H (h);

    Eigen::MatrixXd T (3, 3);
    T << 1,  0,  0,
         0, -2,  0,
         0,  0,  3;
    Eigen::MatrixXd T_inverse = T.inverse();


    H.transform(T);
    H.transform(T_inverse);

    BOOST_CHECK(H.get_matrix_representation().isApprox(h, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_rotate_throws ) {

    // Create a random OneElectronOperator
    size_t dim = 3;
    GQCP::OneElectronOperator M (Eigen::MatrixXd::Random(dim, dim));


    // Check if a non-unitary matrix as transformation matrix causes a throw
    Eigen::MatrixXd U (Eigen::MatrixXd::Random(dim, dim));
    BOOST_CHECK_THROW(M.rotate(U), std::invalid_argument);


    // Check if a unitary matrix as transformation matrix is accepted
    M.rotate(Eigen::MatrixXd::Identity(dim, dim));
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_rotate_JacobiRotationParameters ) {

    // Create a random OneElectronOperator
    size_t dim = 5;
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(dim, dim);
    GQCP::OneElectronOperator M1 (m);
    GQCP::OneElectronOperator M2 (m);


    // Check that using a Jacobi transformation (rotation) matrix as U is equal to the custom transformation (rotation)
    // with custom JacobiRotationParameters
    GQCP::JacobiRotationParameters jacobi_rotation_parameters (4, 2, 56.81);

    auto U = GQCP::jacobiRotationMatrix(jacobi_rotation_parameters, dim);


    M1.rotate(jacobi_rotation_parameters);
    M2.rotate(U);


    BOOST_CHECK(M1.get_matrix_representation().isApprox(M2.get_matrix_representation(), 1.0e-12));
}
