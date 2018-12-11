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
#define BOOST_TEST_MODULE "HoppingMatrix"


#include "HoppingMatrix.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( constructor ) {

    Eigen::MatrixXd H1 (3, 4);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix H (H1), std::invalid_argument);  // not square

    Eigen::MatrixXd H2 = Eigen::MatrixXd::Random(3, 3);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix H (H2), std::invalid_argument);  // not symmetric


    Eigen::MatrixXd H3 = H2 + H2.transpose();
    BOOST_CHECK_NO_THROW(GQCP::HoppingMatrix H (H3));
}


BOOST_AUTO_TEST_CASE ( FromUpperTriangle_throw ) {

    Eigen::VectorXd v1 (8);
    BOOST_CHECK_THROW(GQCP::HoppingMatrix H (v1), std::invalid_argument);


    Eigen::VectorXd v2 (6);
    GQCP::HoppingMatrix H (v2);
    BOOST_CHECK_NO_THROW(GQCP::HoppingMatrix H (v2));
}
