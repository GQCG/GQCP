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
#define BOOST_TEST_MODULE "Eigenpair"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


#include "optimization/Eigenpair.hpp"



BOOST_AUTO_TEST_CASE ( isEqual ) {

    // Create some test eigenvectors
    Eigen::VectorXd u1 (2);
    u1 << 1, 0;
    Eigen::VectorXd u2 (2);
    u2 << -1, 0;
    Eigen::VectorXd v (2);
    v << 1, 1;
    Eigen::VectorXd w (3);
    w << 1, 0, 0;


    // Test isEqual for some eigenpairs
    GQCP::Eigenpair eigenpair1 (1, u1);
    GQCP::Eigenpair eigenpair2 (1, u2);
    GQCP::Eigenpair eigenpair3 (2, u1);
    GQCP::Eigenpair eigenpair4 (1, v);
    GQCP::Eigenpair eigenpair5 (1, w);

    BOOST_CHECK(eigenpair1.isEqual(eigenpair1));
    BOOST_CHECK(eigenpair1.isEqual(eigenpair2));  // eigenvectors are equal up to their sign

    BOOST_CHECK(!eigenpair1.isEqual(eigenpair3));  // 1 != 2
    BOOST_CHECK(!eigenpair1.isEqual(eigenpair4));  // u1 != v


    BOOST_CHECK_THROW(eigenpair1.isEqual(eigenpair5), std::invalid_argument);  // can't compare eigenvectors of different dimensions
}
