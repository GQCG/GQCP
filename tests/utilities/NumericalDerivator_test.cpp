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
#define BOOST_TEST_MODULE "NumericalDerivator"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "utilities/NumericalDerivativeCalculator.hpp"



BOOST_AUTO_TEST_CASE ( derive_xcubed ) {

    GQCP::UnaryFunction xcubed = [](double x) { return pow(x, 3);};
    GQCP::NumericalDerivativeCalculator<4> derivator (xcubed, 0, 0.001);

    BOOST_CHECK(derivator.get_derivative(0) == 0);
    BOOST_CHECK(std::abs(derivator.get_derivative(1) - 1e-06) < 1e-10);
    BOOST_CHECK(std::abs(derivator.get_derivative(2) - 0.006) < 1e-10);
    BOOST_CHECK(std::abs(derivator.get_derivative(3) - 6) < 1e-10);
    BOOST_CHECK(std::abs(derivator.get_derivative(4) - 0) < 1e-10);
}



BOOST_AUTO_TEST_CASE ( derive_eigenproblem ) {

    GQCP::UnaryFunction xcubed = [](double x) { return pow(x, 3);};
    GQCP::NumericEigenProblem example = [xcubed](double x, const Eigen::VectorXd& q) {return GQCP::Eigenpair(xcubed(x),2*q);};
    Eigen::VectorXd q = Eigen::VectorXd::Ones(1);
    GQCP::NumericalGuessDerivativeCalculator<4> derivator (example, 0, 0.001, q);

    BOOST_CHECK(derivator.get_derivative(0) == 0);
    BOOST_CHECK(std::abs(derivator.get_derivative(1) - 1e-06) < 1e-10);
    BOOST_CHECK(std::abs(derivator.get_derivative(2) - 0.006) < 1e-10);
    BOOST_CHECK(std::abs(derivator.get_derivative(3) - 6) < 1e-10);
    BOOST_CHECK(std::abs(derivator.get_derivative(4) - 0) < 1e-10);

    BOOST_CHECK(derivator.get_eigenvector(0).isApprox(2*q));
    BOOST_CHECK(derivator.get_eigenvector(1).isApprox(4*q));
    BOOST_CHECK(derivator.get_eigenvector(2).isApprox(8*q));
    BOOST_CHECK(derivator.get_eigenvector(3).isApprox(16*q));
    BOOST_CHECK(derivator.get_eigenvector(4).isApprox(32*q));
}