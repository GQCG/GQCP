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
#ifndef GQCP_NEWTONSYSTEMOFEQUATIONSSOLVER_HPP
#define GQCP_NEWTONSYSTEMOFEQUATIONSSOLVER_HPP

#include <Eigen/Dense>

#include "optimization/BaseSystemOfEquationsSolver.hpp"
#include "common.hpp"



namespace GQCP {


/**
 *  A class that solves a system of equations by using a Newton algorithm
 */
class NewtonSystemOfEquationsSolver : public BaseSystemOfEquationsSolver {
private:
    VectorFunction f;  // function wrapper for the vector 'function'
    MatrixFunction J;  // function wrapper for the JacobianFunction


public:
    // CONSTRUCTORS
    /**
     *  @param x0                               the initial guess
     *  @param f                                a callable vector function
     *  @param J                                the corresponding callable Jacobian
     *  @param convergence_threshold            the threshold used to determine convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations in the algorithm
     */
    NewtonSystemOfEquationsSolver(const Eigen::VectorXd& x0, const VectorFunction& f, const MatrixFunction& J, double convergence_threshold = 1.0e-08, size_t maximum_number_of_iterations = 128);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  Find a solution to the problem f(x) = 0
     *
     *  If successful, it sets
     *      - is_solved to true
     *      - the found solution
     */
    void solve() override;
};


}  // namespace GQCP


#endif  // GQCP_NEWTONSYSTEMOFEQUATIONSSOLVER_HPP
