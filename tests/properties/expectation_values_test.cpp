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
#define BOOST_TEST_MODULE "expectation_values"

#include "properties/expectation_values.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( one_electron_throw ) {

    GQCP::OneElectronOperator h (Eigen::MatrixXd::Zero(2, 2));
    GQCP::OneRDM D_valid (Eigen::MatrixXd::Zero(2, 2));
    GQCP::OneRDM D_invalid (Eigen::MatrixXd::Zero(3, 3));

    BOOST_CHECK_THROW(GQCP::calculateExpectationValue(h, D_invalid), std::invalid_argument);
    BOOST_CHECK_NO_THROW(GQCP::calculateExpectationValue(h, D_valid));
}


BOOST_AUTO_TEST_CASE ( two_electron_throw ) {

    Eigen::Tensor<double, 4> g_tensor (2, 2, 2, 2);
    g_tensor.setZero();
    GQCP::TwoElectronOperator g (g_tensor);

    Eigen::Tensor<double, 4> d_tensor_valid (2, 2, 2, 2);
    d_tensor_valid.setZero();
    Eigen::Tensor<double, 4> d_tensor_invalid (3, 3, 3, 3);
    d_tensor_valid.setZero();
    GQCP::TwoRDM d_valid (d_tensor_valid);
    GQCP::TwoRDM d_invalid (d_tensor_invalid);

    BOOST_CHECK_THROW(GQCP::calculateExpectationValue(g, d_invalid), std::invalid_argument);
    BOOST_CHECK_NO_THROW(GQCP::calculateExpectationValue(g, d_valid));
}
