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
#define BOOST_TEST_MODULE "SquareFourIndexTensor"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "math/SquareFourIndexTensor.hpp"


BOOST_AUTO_TEST_CASE ( constructor ) {

    Eigen::Tensor<double, 4> T1 (2, 2, 2, 2);
    T1.setZero();

    BOOST_CHECK_NO_THROW(GQCP::SquareFourIndexTensor<double> square_T1 (T1));

    Eigen::Tensor<double, 4> T2 (2, 1, 2, 2);
    BOOST_CHECK_THROW(GQCP::SquareFourIndexTensor<double> square_T2 (T1), std::invalid_argument);
}
