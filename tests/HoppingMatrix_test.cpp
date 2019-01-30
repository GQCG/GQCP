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
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "HoppingMatrix.hpp"


BOOST_AUTO_TEST_CASE ( constructor_throws ) {

    Eigen::MatrixXd H1 (3, 4);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix H (H1), std::invalid_argument);  // not square

    Eigen::MatrixXd H2 = Eigen::MatrixXd::Random(3, 3);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix H (H2), std::invalid_argument);  // not symmetric


    Eigen::MatrixXd H3 = H2 + H2.transpose();
    BOOST_CHECK_NO_THROW(GQCP::HoppingMatrix H (H3));
}


BOOST_AUTO_TEST_CASE ( FromUpperTriangle_throws ) {

    Eigen::VectorXd v1 (8);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix::FromUpperTriangle(v1), std::invalid_argument);  // 8 is not a square number


    Eigen::VectorXd v2 (6);
    BOOST_CHECK_NO_THROW(GQCP::HoppingMatrix::FromUpperTriangle(v2));
}


BOOST_AUTO_TEST_CASE ( triangle_adjacency_matrix ) {

    // Check the conversion to a triangle hopping matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);
    A << 0, 1, 1,
         1, 0, 1,
         1, 1, 0;
    double t = 1.0;
    double U = 2.0;
    GQCP::HoppingMatrix H (A, t, U);


    Eigen::MatrixXd H_ref = Eigen::MatrixXd::Zero(3, 3);
    H_ref << U, -t, -t,
            -t,  U, -t,
            -t, -t,  U;


    BOOST_CHECK(H_ref.isApprox(H.asMatrix()));
}
