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
#define BOOST_TEST_MODULE "ChemicalMatrix"

#include <boost/test/unit_test.hpp>

#include "Mathematical/ChemicalMatrix.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check the constructor API for ChemicalMatrix
 */
BOOST_AUTO_TEST_CASE ( constructor ) {

    const GQCP::ChemicalMatrix<double> M1 = GQCP::ChemicalMatrix<double>::Zero(2, 2);
    BOOST_CHECK_NO_THROW(GQCP::ChemicalMatrix<double> chemical_matrix (M1));

    const GQCP::MatrixX<double> M2 = GQCP::MatrixX<double>::Zero(2, 1);
    BOOST_CHECK_THROW(GQCP::ChemicalMatrix<double> chemical_matrix (M2), std::invalid_argument);
}


/**
 *  Check the basis transformation formula with a trivial transformation: T being a unit matrix
 */
BOOST_AUTO_TEST_CASE ( ChemicalMatrix_transform_trivial ) {

    const GQCP::ChemicalMatrix<double> T = GQCP::ChemicalMatrix<double>::Identity(3, 3);

    GQCP::ChemicalMatrix<double> h = GQCP::ChemicalMatrix<double>::Random(3, 3);
    const GQCP::ChemicalMatrix<double> h_copy = h;
    h.basisTransformInPlace(T);

    BOOST_CHECK(h_copy.isApprox(h, 1.0e-12));
}


/**
 *  Check if we transform with a transformation matrix and its inverse, we get a zero operation
 */
BOOST_AUTO_TEST_CASE ( ChemicalMatrix_transform_and_inverse ) {

    GQCP::ChemicalMatrix<double> T (3);
    T << 1,  0,  0,
         0, -2,  0,
         0,  0,  3;
    const GQCP::ChemicalMatrix<double> T_inverse = T.inverse();


    GQCP::ChemicalMatrix<double> h = GQCP::ChemicalMatrix<double>::Random(3, 3);
    GQCP::ChemicalMatrix<double> h_copy = h;
    h.basisTransformInPlace(T);
    h.basisTransformInPlace(T_inverse);

    BOOST_CHECK(h.isApprox(h_copy, 1.0e-12));
}


/**
 *  Check if we can't rotate with a unitary matrix
 */
BOOST_AUTO_TEST_CASE ( ChemicalMatrix_rotate_throws ) {

    // Create a random ChemicalMatrix
    size_t dim = 3;
    GQCP::ChemicalMatrix<double> M = GQCP::ChemicalMatrix<double>::Random(dim, dim);


    // Check if a non-unitary matrix as transformation matrix causes a throw
    const GQCP::SquareMatrix<double> T = GQCP::SquareMatrix<double>::Random(dim, dim);
    BOOST_CHECK_THROW(M.basisRotateInPlace(T), std::invalid_argument);


    // Check if a unitary matrix as transformation matrix is accepted
    const GQCP::SquareMatrix<double> U= GQCP::SquareMatrix<double>::Identity(dim, dim);
    M.basisRotateInPlace(U);
}


/**
 *  Check if rotating with JacobiRotationParameters is the same as with the corresponding Jacobi rotation matrix
 */
BOOST_AUTO_TEST_CASE ( SQOneElectronOperator_rotate_JacobiRotationParameters ) {

    // Create a random ChemicalMatrix
    const size_t dim = 5;
    GQCP::ChemicalMatrix<double> M1 = GQCP::ChemicalMatrix<double>::Random(dim, dim);
    auto M2 = M1;

    // Create random Jacobi rotation parameters and the corresponding Jacobi rotation matrix
    GQCP::JacobiRotationParameters jacobi_rotation_parameters (4, 2, 56.81);
    const auto J = GQCP::SquareMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);


    // Rotate and check the result
    M1.basisRotateInPlace(jacobi_rotation_parameters);
    M2.basisRotateInPlace(J);

    BOOST_CHECK(M1.isApprox(M2, 1.0e-12));
}
