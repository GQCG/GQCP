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

#define BOOST_TEST_MODULE "QCMatrix"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/QCMatrix.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check the constructor API for QCMatrix
 */
BOOST_AUTO_TEST_CASE(constructor) {

    const GQCP::QCMatrix<double> M1 = GQCP::QCMatrix<double>::Zero(2, 2);
    BOOST_CHECK_NO_THROW(GQCP::QCMatrix<double> chemical_matrix {M1});

    const GQCP::MatrixX<double> M2 = GQCP::MatrixX<double>::Zero(2, 1);
    BOOST_CHECK_THROW(GQCP::QCMatrix<double> chemical_matrix {M2}, std::invalid_argument);
}


/**
 *  Check the basis transformation formula with a trivial transformation: T being a unit matrix
 */
BOOST_AUTO_TEST_CASE(QCMatrix_transform_trivial) {

    const GQCP::QCMatrix<double> T = GQCP::QCMatrix<double>::Identity(3, 3);

    GQCP::QCMatrix<double> h = GQCP::QCMatrix<double>::Random(3, 3);
    const GQCP::QCMatrix<double> h_copy = h;
    h.basisTransform(T);

    BOOST_CHECK(h_copy.isApprox(h, 1.0e-12));
}


/**
 *  Check if we transform with a transformation matrix and its inverse, we get a zero operation
 */
BOOST_AUTO_TEST_CASE(QCMatrix_transform_and_inverse) {

    GQCP::QCMatrix<double> T {3};
    // clang-format off
    T << 1,  0,  0,
         0, -2,  0,
         0,  0,  3;
    // clang-format on
    const GQCP::QCMatrix<double> T_inverse = T.inverse();


    GQCP::QCMatrix<double> h = GQCP::QCMatrix<double>::Random(3, 3);
    GQCP::QCMatrix<double> h_copy = h;
    h.basisTransform(T);
    h.basisTransform(T_inverse);

    BOOST_CHECK(h.isApprox(h_copy, 1.0e-12));
}


/**
 *  Check if we can't rotate with a unitary matrix
 */
BOOST_AUTO_TEST_CASE(QCMatrix_rotate_throws) {

    // Create a random QCMatrix
    size_t dim = 3;
    GQCP::QCMatrix<double> M = GQCP::QCMatrix<double>::Random(dim, dim);


    // Check if a non-unitary matrix as transformation matrix causes a throw
    const GQCP::TransformationMatrix<double> T = GQCP::TransformationMatrix<double>::Random(dim);
    BOOST_CHECK_THROW(M.basisRotate(T), std::invalid_argument);


    // Check if a unitary matrix as transformation matrix is accepted
    const GQCP::TransformationMatrix<double> U = GQCP::TransformationMatrix<double>::Identity(dim);
    M.basisRotate(U);
}


/**
 *  Check if rotating with JacobiRotationParameters is the same as with the corresponding Jacobi rotation matrix
 */
BOOST_AUTO_TEST_CASE(SQOneElectronOperator_rotate_JacobiRotationParameters) {

    // Create a random QCMatrix
    const size_t dim = 5;
    GQCP::QCMatrix<double> M1 = GQCP::QCMatrix<double>::Random(dim, dim);
    auto M2 = M1;

    // Create random Jacobi rotation parameters and the corresponding Jacobi rotation matrix
    GQCP::JacobiRotationParameters jacobi_rotation_parameters {4, 2, 56.81};
    const auto J = GQCP::TransformationMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);


    // Rotate and check the result
    M1.basisRotate(jacobi_rotation_parameters);
    M2.basisRotate(J);

    BOOST_CHECK(M1.isApprox(M2, 1.0e-12));
}


/**
 *  Check if a unitary transformation (i.e. a rotation) leaves the overlap matrix invariant
 */
BOOST_AUTO_TEST_CASE(rotate_overlap_invariant) {

    // Initialize the overlap matrix
    const size_t K = 3;
    GQCP::QCMatrix<double> S {K};
    // clang-format off
    S << 1.0, 0.5, 0.0,
         0.5, 2.0, 0.0,
         0.0, 0.0, 1.0;
    // clang-format on


    // Initialize the reference
    GQCP::QCMatrix<double> S_rotated_ref {K};
    // clang-format off
    S_rotated_ref <<  1.0, 0.0, -0.5,
                      0.0, 1.0,  0.0,
                     -0.5, 0.0,  2.0;
    // clang-format on


    // Rotate the overlap matrix and check the result
    GQCP::TransformationMatrix<double> U {K};
    // clang-format off
    U << 1.0,  0.0,  0.0,
         0.0,  0.0, -1.0,
         0.0, -1.0,  0.0;
    // clang-format on

    auto S_copy = S;
    S_copy.basisRotate(U);
    BOOST_CHECK(S_copy.isApprox(S_rotated_ref, 1.0e-08));
}


/**
 *  Check a rotation through Jacobi rotation parameters through a manual calculation
 */
BOOST_AUTO_TEST_CASE(rotate_Jacobi_manual) {

    // Initialize the test overlap matrix and rotation parameters
    const size_t K = 3;
    GQCP::QCMatrix<double> S {K};
    // clang-format off
    S << 1.0, 0.5, 0.0,
         0.5, 2.0, 0.0,
         0.0, 0.0, 1.0;
    // clang-format on

    const GQCP::JacobiRotationParameters jacobi_rotation_parameters {1, 0, boost::math::constants::half_pi<double>()};  // interchanges two orbitals


    // Initialize the reference, rotate the overlap matrix and check the result
    GQCP::QCMatrix<double> S_rotated_ref {K};  // manual calculation
    // clang-format off
    S_rotated_ref <<  2.0, -0.5, 0.0,
                     -0.5,  1.0, 0.0,
                      0.0,  0.0, 1.0;
    // clang-format on


    S.basisRotate(jacobi_rotation_parameters);
    BOOST_CHECK(S.isApprox(S_rotated_ref, 1.0e-08));
}
