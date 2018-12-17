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
#define BOOST_TEST_MODULE "OneRDM"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "RDM/OneRDM.hpp"


BOOST_AUTO_TEST_CASE ( OneRDM_constructor ) {

    // Check a correct constructor
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(4, 4);
    GQCP::OneRDM D (matrix);


    // Check a faulty constructor
    Eigen::MatrixXd matrix2 = Eigen::MatrixXd::Zero(3, 4);
    BOOST_CHECK_THROW(GQCP::OneRDM D2 (matrix2), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE ( isEqualTo ) {

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3, 3);
    GQCP::OneRDM O1 (A);
    GQCP::OneRDM O2 (A);
    BOOST_CHECK(O1.isEqualTo(O2, 1.0e-05));

    GQCP::OneRDM O3 (2*A);
    BOOST_CHECK(!(O3.isEqualTo(O1)));
}


BOOST_AUTO_TEST_CASE ( operator_equals ) {

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3, 3);
    GQCP::OneRDM O1 (A);
    GQCP::OneRDM O2 (A);
    BOOST_CHECK(O1 == O2);
}
