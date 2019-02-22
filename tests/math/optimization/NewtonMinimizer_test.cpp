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
#define BOOST_TEST_MODULE "NewtonMinimizer"
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


#include "math/optimization/NewtonMinimizer.hpp"



/*
 *  TEST FUNCTION DEFINITIONS
 */

/**
 *  Implement a simple scalar function that returns x.x
 */
double f(const Eigen::VectorXd& x) {
    return x.squaredNorm();
}


/**
 *  Implement the gradient of the scalar function
 */
Eigen::VectorXd grad(const Eigen::VectorXd& x) {
    return 2 * x;
}


/**
 *  Implement the Hessian of the scalar function
 */
Eigen::MatrixXd H(const Eigen::VectorXd& x) {
    return 2 * Eigen::MatrixXd::Identity(x.size(), x.size());
}



/*
 *  BOOST UNIT TESTS
 */

BOOST_AUTO_TEST_CASE ( norm_squared_function_minimization ) {

    // Test that the previous implementations actually work by checking the values at x=(1,1)
    Eigen::VectorXd x_test (2);
    x_test << 1, 1;

    Eigen::VectorXd grad_test (2);
    grad_test << 2, 2;

    Eigen::MatrixXd H_test (2, 2);
    H_test << 2, 0,
              0, 2;

    BOOST_REQUIRE(std::abs(f(x_test) - 2.0) < 1.0e-12);
    BOOST_REQUIRE(grad_test.isApprox(grad(x_test), 1.0e-8));
    BOOST_REQUIRE(H_test.isApprox(H(x_test), 1.0e-8));


    // Do the numerical optimization
    Eigen::VectorXd x0 (2);
    x0 << 4, 2;

    GQCP::NewtonMinimizer newton_minimizer (x0, grad, H);  // apparently, the compiler can convert to numopt::VectorFunction and numopt::JacobianFunction
    newton_minimizer.solve();
    Eigen::VectorXd solution = newton_minimizer.get_solution();


    BOOST_CHECK(solution.isZero(1.0e-08));  // the analytical minimizer of f(x) is x=(0,0)
}
