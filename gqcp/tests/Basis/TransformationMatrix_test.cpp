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

#define BOOST_TEST_MODULE "TransformationMatrix_test"

#include <boost/test/unit_test.hpp>

#include "Basis/TransformationMatrix.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check if the 'product' of two transformation matrices behaves as expected
 */
BOOST_AUTO_TEST_CASE(transform) {

    // Create two transformation matrices
    GQCP::TransformationMatrix<double> T1 {2};
    // clang-format off
    T1 << 1.0, 0.0,
          1.0, 3.0;
    // clang-format on

    GQCP::TransformationMatrix<double> T2 {2};
    // clang-format off
    T2 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on


    // Set up and check the expected result
    GQCP::TransformationMatrix<double> T_total {2};
    // clang-format off
    T_total <<  1.0,  2.0,
               10.0, 14.0;
    // clang-format on

    T1.transform(T2);  // first do T1, then T2
    BOOST_CHECK(T1.isApprox(T_total, 1.0e-08));
}


/**
 *  Check the construction of Jacobi rotation matrices from JacobiRotationParameters.
 */
BOOST_AUTO_TEST_CASE(FromJacobi) {

    // A random Jacobi matrix is unitary
    BOOST_CHECK(GQCP::TransformationMatrix<double>::FromJacobi(GQCP::JacobiRotationParameters(7, 4, 6.9921), 10).isUnitary());
    BOOST_CHECK(GQCP::TransformationMatrix<double>::FromJacobi(GQCP::JacobiRotationParameters(9, 1, 78.00166), 22).isUnitary());


    // Check the cos, sin, -sin, cos convention. Since p>q, this means we should end up with cos, -sin, sin, cos in the Jacobi rotation matrix
    const auto half_pi = boost::math::constants::half_pi<double>();
    const auto J = GQCP::TransformationMatrix<double>::FromJacobi(GQCP::JacobiRotationParameters(1, 0, half_pi), 2);

    BOOST_CHECK(std::abs(J(0, 0) - 0) < 1.0e-12);
    BOOST_CHECK(std::abs(J(0, 1) - (-1)) < 1.0e-12);
    BOOST_CHECK(std::abs(J(1, 0) - 1) < 1.0e-12);
    BOOST_CHECK(std::abs(J(0, 0) - 0) < 1.0e-12);
}
