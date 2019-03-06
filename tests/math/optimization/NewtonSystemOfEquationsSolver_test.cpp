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
#define BOOST_TEST_MODULE "NewtonSystemOfEquations"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "math/optimization/NewtonSystemOfEquationsSolver.hpp"
#include "math/SquareMatrix.hpp"


/*
 *  TEST FUNCTION DEFINITIONS
 */

/**
 *  Implement a simple vector function that returns (x.x, 2x.x)
 */
GQCP::VectorX<double> f(const GQCP::VectorX<double>& x) {

    GQCP::VectorX<double> f (2);

    f << x(0) * x(0) + x(1) * x(1),
         2 * x(0) * x(0) + 2 * x(1) * x(1);

    return f;
}


/**
 *  Implement the Jacobian of the previous function
 */
GQCP::SquareMatrix<double> J(const GQCP::VectorX<double>& x) {

    GQCP::SquareMatrix<double> J (2);

    J << 2 * x(0), 2 * x(1),
         4 * x(0), 4 * x(1);

    return J;
}



/*
 *  BOOST UNIT TESTS
 */

BOOST_AUTO_TEST_CASE ( norm_squared_function_syseq ) {

    // Test that the previous implementations actually work by checking the values at x=(1,1)
    GQCP::VectorX<double> f_test (2);
    f_test << 2, 4;

    GQCP::SquareMatrix<double> J_test (2);
    J_test << 2, 2,
              4, 4;

    GQCP::VectorX<double> x_test (2);
    x_test << 1, 1;

    BOOST_REQUIRE(f_test.isApprox(f(x_test), 1.0e-8));
    BOOST_REQUIRE(J_test.isApprox(J(x_test), 1.0e-8));


    // Do the numerical optimization
    GQCP::VectorX<double> x0 (2);
    x0 << 3, 2;
    GQCP::NewtonSystemOfEquationsSolver newton_vector_opt (x0, f, J);  // apparently, the compiler can convert to GQCP::VectorFunction and GQCP::JacobianFunction
    newton_vector_opt.solve();
    GQCP::VectorX<double> solution = newton_vector_opt.get_solution();


    BOOST_CHECK(solution.isZero(1.0e-08));  // the analytical solution of f(x) = (0,0) is x=(0,0)
}
