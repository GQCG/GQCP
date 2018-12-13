// This file is part of GQCG-numopt.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-numopt is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-numopt is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-numopt.  If not, see <http://www.gnu.org/licenses/>.
#include "NewtonMinimizer.hpp"

#include "common.hpp"
#include "NewtonSystemOfEquationsSolver.hpp"



#include <iostream>


namespace numopt {
namespace minimization {


/*
 *  CONSTRUCTORS
 */
/**
 *  @param x0                           the initial guess
 *  @param grad                         the callable gradient function
 *  @param H                            the callable Hessian function
 *  @param convergence_threshold        the threshold used for establishing convergence
 */
NewtonMinimizer::NewtonMinimizer(const Eigen::VectorXd& x0, const VectorFunction& grad, const MatrixFunction& H,
                                 double convergence_threshold) :
    BaseMinimizer(x0, convergence_threshold),
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

    // Requiring the gradient to vanish, and then calculating the corresponding Newton step, is equivalent to solving
    // the system of equations grad(f(x)) = 0 using a Newton step

    // For mathematical correctness, the Jacobian of the gradient is the transpose of the Hessian of the scalar function
    // behind it
    numopt::MatrixFunction H_t = [this](const Eigen::VectorXd& x) {
        Eigen::MatrixXd H = this->H(x);
        H.transposeInPlace();
        return H;
    };
    // We have defined an elaborate lambda function, because
    // numopt::JacobianFunction H_t = [this](const Eigen::VectorXd& x) { return this->H(x).transpose(); }
    // produces an error


    // For previously established reasons, we can use the NewtonSystemOfEquationsSolver as an implementation of this
    // minimization problem
    numopt::syseq::NewtonSystemOfEquationsSolver newton_syseq_solver (this->x0, this->grad, H_t, this->convergence_threshold);
    newton_syseq_solver.solve();


    // If we haven't found a solution, the error is raised inside the NewtonSystemOfEquationsSolver, so we are free to
    // assert that if the data flow gets to here, a solution is found
    this->is_solved = true;
    this->x = newton_syseq_solver.get_solution();
}


}  // namespace minimization
}  // namespace numopt
