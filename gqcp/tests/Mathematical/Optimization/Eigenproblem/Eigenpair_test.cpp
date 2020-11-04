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

#define BOOST_TEST_MODULE "Eigenpair"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Eigenproblem/Eigenpair.hpp"


BOOST_AUTO_TEST_CASE(isEqualTo) {

    // Create some test eigenvectors
    GQCP::VectorX<double> u1 {2};
    u1 << 1, 0;
    GQCP::VectorX<double> u2 {2};
    u2 << -1, 0;
    GQCP::VectorX<double> v {2};
    v << 1, 1;
    GQCP::VectorX<double> w {3};
    w << 1, 0, 0;


    // Test isEqualTo for some eigenpairs
    GQCP::Eigenpair<double> eigenpair1 {1, u1};
    GQCP::Eigenpair<double> eigenpair2 {1, u2};
    GQCP::Eigenpair<double> eigenpair3 {2, u1};
    GQCP::Eigenpair<double> eigenpair4 {1, v};
    GQCP::Eigenpair<double> eigenpair5 {1, w};

    BOOST_CHECK(eigenpair1.isEqualTo(eigenpair1));
    BOOST_CHECK(eigenpair1.isEqualTo(eigenpair2));  // eigenvectors are equal up to their sign

    BOOST_CHECK(!eigenpair1.isEqualTo(eigenpair3));  // 1 != 2
    BOOST_CHECK(!eigenpair1.isEqualTo(eigenpair4));  // u1 != v


    BOOST_CHECK_THROW(eigenpair1.isEqualTo(eigenpair5), std::invalid_argument);  // can't compare eigenvectors of different dimensions
}
