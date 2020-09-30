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

#define BOOST_TEST_MODULE "SimpleTransformationMatrix"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/RTransformationMatrix.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check if the Jacobi rotation is implemented correctly.
 */
BOOST_AUTO_TEST_CASE(jacobi_rotation) {

    const size_t dim = 3;

    // Set up an identity transformation matrix.
    GQCP::RTransformationMatrix<double> T = GQCP::SquareMatrix<double>::Identity(dim);
    const GQCP::JacobiRotationParameters jacobi {1, 0, boost::math::constants::half_pi<double>()};


    // Set up the reference and check the results.
    GQCP::SquareMatrix<double> T_rotated_ref {dim};
    // clang-format off
    T_rotated_ref << 0.0, -1.0, 0.0,
                     1.0,  0.0, 0.0,
                     0.0,  0.0, 1.0;
    // clang-format on


    // Check the in-place and the returning methods.
    const auto T_transformed = T.rotated(jacobi);
    BOOST_CHECK(T_transformed.isApprox(T_rotated_ref, 1.0e-12));

    T.rotate(jacobi);
    BOOST_CHECK(T.isApprox(T_rotated_ref, 1.0e-12));
}


/**
 *  Check if Jacobi transformation matrices are unitary.
 */
BOOST_AUTO_TEST_CASE(FromJacobi_unitary) {

    BOOST_CHECK(GQCP::RTransformationMatrix<double>::FromJacobi(GQCP::JacobiRotationParameters(7, 4, 6.9921), 10).isUnitary());
    BOOST_CHECK(GQCP::RTransformationMatrix<double>::FromJacobi(GQCP::JacobiRotationParameters(9, 1, 78.00166), 22).isUnitary());
}


/**
 *  Check if the 'product' of two transformation matrices behaves as expected.
 */
BOOST_AUTO_TEST_CASE(transform) {

    // Create two transformation matrices.
    GQCP::RTransformationMatrix<double> T1 {2};
    // clang-format off
    T1 << 1.0, 0.0,
          1.0, 3.0;
    // clang-format on

    GQCP::RTransformationMatrix<double> T2 {2};
    // clang-format off
    T2 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on


    // Set up and check the expected result
    GQCP::RTransformationMatrix<double> T_total_ref {2};
    // clang-format off
    T_total_ref <<  1.0,  2.0,
                   10.0, 14.0;
    // clang-format on


    // Check the in-place and the returning methods.
    const auto T_total = T1.transformed(T2);
    BOOST_CHECK(T_total.isApprox(T_total_ref, 1.0e-12));

    T1.transform(T2);
    BOOST_CHECK(T1.isApprox(T_total_ref, 1.0e-12));
}
