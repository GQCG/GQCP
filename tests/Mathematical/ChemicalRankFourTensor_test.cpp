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
#define BOOST_TEST_MODULE "ChemicalRankFourTensor"

#include <boost/test/unit_test.hpp>

#include "Mathematical/ChemicalRankFourTensor.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check the interface for ChemicalRankFourTensor constructors
 */
BOOST_AUTO_TEST_CASE ( constructor ) {

    GQCP::Tensor<double, 4> T1 (2, 2, 2, 2);
    T1.setZero();
    BOOST_CHECK_NO_THROW(GQCP::ChemicalRankFourTensor<double> square_T1 (T1));

    const GQCP::Tensor<double, 4> T2 (2, 1, 2, 2);
    BOOST_CHECK_THROW(GQCP::ChemicalRankFourTensor<double> square_T2 (T2), std::invalid_argument);  // not square
}


/**
 *  Check the basis transformation formula for a trivial case: T being a unit matrix
 */
BOOST_AUTO_TEST_CASE ( ChemicalRankFourTensor_basisTransformInPlace_trivial ) {

    const GQCP::SquareMatrix<double> T = GQCP::SquareMatrix<double>::Identity(3, 3);

    GQCP::ChemicalRankFourTensor<double> G (3);
    const auto G_copy = G;  // the reference
    G.basisTransformInPlace(T);

    BOOST_CHECK(G_copy.isApprox(G, 1.0e-12));
}


/**
 *  Check the basis transformation formula using an other implementation (the old olsens code) from Ayers' Lab
 */
BOOST_AUTO_TEST_CASE ( SQTwoElectronOperator_transform_olsens ) {

    // Set an example transformation matrix and two-electron integrals tensor
    const size_t dim = 2;
    GQCP::SquareMatrix<double> T (dim);
    T << 1, 2,
         3, 4;

    GQCP::ChemicalRankFourTensor<double> G (dim);
    G.setZero();
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    G(i, j, k, l) = l + 2*k + 4*j + 8*i;
                }
            }
        }
    }
    G.basisTransformInPlace(T);


    // Read in the reference and check
    const GQCP::ChemicalRankFourTensor<double> g_transformed_ref = GQCP::ChemicalRankFourTensor<double>::FromFile("data/rotated_two_electron_integrals_olsens.data", dim);
    BOOST_CHECK(G.isApprox(g_transformed_ref, 1.0e-12));
}


/**
 *  Check if the code throws errors concerning non-unitary matrices being given to .rotate()
 */
BOOST_AUTO_TEST_CASE ( ChemicalRankFourTensor_rotate_throws ) {

    // Create a random ChemicalRankFourTensor
    const size_t dim = 3;
    GQCP::ChemicalRankFourTensor<double> g (dim);
    g.setRandom();


    // Check if a non-unitary matrix as transformation matrix causes a throw
    const GQCP::SquareMatrix<double> T = GQCP::SquareMatrix<double>::Random(dim, dim);  // chances are practically zero that a random matrix is unitary
    BOOST_CHECK_THROW(g.basisRotateInPlace(GQCP::SquareMatrix<double>(T)), std::invalid_argument);


    // Check if a unitary matrix as transformation matrix is accepted
    GQCP::SquareMatrix<double> U = GQCP::SquareMatrix<double>::Identity(dim, dim);
    g.basisRotateInPlace(U);
}


/**
 *  Check that rotating using JacobiRotationParameters is the same as using the corresponding Jacobi rotation matrix
 */
BOOST_AUTO_TEST_CASE ( ChemicalRankFourTensor_basisRotateInPlace_JacobiRotationParameters ) {

    // Create a random ChemicalRankFourTensor
    const size_t dim = 5;
    GQCP::ChemicalRankFourTensor<double> g (dim);
    g.setRandom();
    GQCP::ChemicalRankFourTensor<double> G1 (g);
    GQCP::ChemicalRankFourTensor<double> G2 (g);


    const GQCP::JacobiRotationParameters jacobi_rotation_parameters (4, 2, 56.81);
    const auto U = GQCP::SquareMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);

    G1.basisRotateInPlace(jacobi_rotation_parameters);
    G2.basisRotateInPlace(U);


    BOOST_CHECK(G1.isApprox(G2, 1.0e-12));
}
