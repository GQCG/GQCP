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
#include "Mathematical/Optimization/Minimization/NewtonMinimizer.hpp"


#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationSolver.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Utilities/typedefs.hpp"

#include <iostream>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */
/**
 *  @param x0                               the initial guess
 *  @param grad                             the callable gradient function
 *  @param H                                the callable Hessian function
 *  @param convergence_threshold            the threshold used for establishing convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations in the algorithm
 */
NewtonMinimizer::NewtonMinimizer(const VectorX<double>& x0, const VectorFunction<double>& grad, const MatrixFunction<double>& H, double convergence_threshold, size_t maximum_number_of_iterations) :
    BaseMinimizer(x0, convergence_threshold, maximum_number_of_iterations),
    grad (grad),
    H (H)
{}



/*
 *  PUBLIC OVERRIDDEN FUNCTIONS
 */
/**
 *  Minimize the function f(x)
 *
 *  If successful, sets
 *      - is_solved to true
 *      - the found solution
 */
void NewtonMinimizer::solve() {

    // Requiring the gradient to vanish, and then calculating the corresponding Newton step, is equivalent to solving the system of equations grad(f(x)) = 0 using a Newton step

    // For mathematical correctness, the Jacobian of the gradient is the transpose of the Hessian of the scalar function behind it
    MatrixFunction<double> H_t = [this](const VectorX<double>& x) {
        SquareMatrix<double> H = this->H(x);
        H.transposeInPlace();
        return H;
    };


    // For previously established reasons, we can use the a Newton-based system of equations solver for this minimization problem
    NonLinearEquationEnvironment<double> non_linear_environment (this->x, this->grad, H_t);
    auto non_linear_solver = NonLinearEquationSolver<double>::Newton(this->convergence_threshold, this->maximum_number_of_iterations);
    non_linear_solver.iterate(non_linear_environment);
    this->x = non_linear_environment.variables.back();

    // If we haven't found a solution, the error is raised inside the iterative algorithm, so we are free to assert that if the data flow gets to here, a solution is found
    this->is_solved = true;
}


}  // namespace GQCP
