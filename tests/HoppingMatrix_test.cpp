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
#define BOOST_TEST_MODULE "HoppingMatrix"


#include <boost/test/unit_test.hpp>

#include "HoppingMatrix.hpp"


BOOST_AUTO_TEST_CASE ( constructor_throws ) {

    GQCP::MatrixX<double> H1 (3, 4);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix H (H1), std::invalid_argument);  // not square

    GQCP::MatrixX<double> H2 = GQCP::MatrixX<double>::Random(3, 3);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix H (H2), std::invalid_argument);  // not symmetric


    GQCP::MatrixX<double> H3 = H2 + H2.transpose();
    BOOST_CHECK_NO_THROW(GQCP::HoppingMatrix H (H3));
}


BOOST_AUTO_TEST_CASE ( FromUpperTriangle_throws ) {

    GQCP::VectorX<double> v1 (8);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix::FromUpperTriangle(v1), std::invalid_argument);  // 8 is not a square number


    GQCP::VectorX<double> v2 (6);
    BOOST_CHECK_NO_THROW(GQCP::HoppingMatrix::FromUpperTriangle(v2));
}


BOOST_AUTO_TEST_CASE ( triangle_adjacency_matrix ) {

    // Check the conversion to a triangle hopping matrix
    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Zero(3, 3);
    A << 0, 1, 1,
         1, 0, 1,
         1, 1, 0;
    double t = 1.0;
    double U = 2.0;
    GQCP::HoppingMatrix H (A, t, U);


    GQCP::SquareMatrix<double> H_ref = GQCP::SquareMatrix<double>::Zero(3, 3);
    H_ref << U, -t, -t,
            -t,  U, -t,
            -t, -t,  U;


    BOOST_CHECK(H_ref.isApprox(H));
}


BOOST_AUTO_TEST_CASE ( FromCSLine ) {

    GQCP::SquareMatrix<double> H_ref = GQCP::SquareMatrix<double>::Zero(3, 3);
    H_ref << 1, 2, 4,
             2, 3, 5,
             4, 5, 6;

    
    const auto H = GQCP::HoppingMatrix::FromCSLine("1,2,4,3,5,6");

    BOOST_CHECK(H_ref.isApprox(H));
}
