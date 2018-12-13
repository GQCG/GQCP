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
#ifndef NUMOPT_NEWTONSYSTEMOFEQUATIONSSOLVER_HPP
#define NUMOPT_NEWTONSYSTEMOFEQUATIONSSOLVER_HPP

#include <Eigen/Dense>

#include "BaseSystemOfEquationsSolver.hpp"
#include "common.hpp"



namespace numopt {
namespace syseq {


/**
 *  A class that is solves a system of equations by using a Newton algorithm
 */
class NewtonSystemOfEquationsSolver : public BaseSystemOfEquationsSolver {
private:
    const VectorFunction f;  // function wrapper for the vector 'function'
    const MatrixFunction J;  // function wrapper for the JacobianFunction


public:
    // CONSTRUCTORS
    /**
     *  @param x0                           the initial guess
     *  @param f                            a callable vector function
     *  @param J                            the corresponding callable Jacobian
     *  @param convergence_threshold        the threshold used to determine convergence
     */
    NewtonSystemOfEquationsSolver(const Eigen::VectorXd& x0, const VectorFunction& f, const MatrixFunction& J,
                                  double convergence_threshold = 1.0e-08);


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


}  // namespace syseq
}  // namespace numopt


#endif // NUMOPT_NEWTONSYSTEMOFEQUATIONSSOLVER_HPP
