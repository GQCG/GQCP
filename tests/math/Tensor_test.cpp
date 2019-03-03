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
#define BOOST_TEST_MODULE "Tensor"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "math/Tensor.hpp"


BOOST_AUTO_TEST_CASE ( constructor_assignment ) {

    // A small check to see if the interface of the constructor and assignment operator works as expected

    Eigen::Tensor<double, 3> A (2, 2, 2);
    Eigen::Tensor<double, 3> B (2, 2, 2);

    GQCP::Tensor<double, 3> T1 (A + B);
    GQCP::Tensor<double, 3> T2 = 2 * B;
}
