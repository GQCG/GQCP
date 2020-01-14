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
#pragma once


#include "Mathematical/Representation/Matrix.hpp"
#include "typedefs.hpp"


namespace GQCP {


/**
 *  A class that solves a system of non-linear equations by using an exact Newton-step algorithm
 */
class NewtonNLSystemOfEquationsSolver {
private:
    double convergence_threshold;  // the threshold on the norm of the update that determines convergence
    size_t maximum_number_of_iterations;

    VectorFunction f;  // a callable vector function that returns the value of each of the equations in a resulting vector, given the current variables
    MatrixFunction J;  // the corresponding callable Jacobian


public:
    // CONSTRUCTORS

    /**
     *  @param f                                a callable vector function that returns the value of each of the equations in a resulting vector, given the current variables
     *  @param J                                the corresponding callable Jacobian
     *  @param convergence_threshold            the threshold on the norm of the update that determines convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations in the algorithm
     */
    NewtonNLSystemOfEquationsSolver(const VectorFunction& f, const MatrixFunction& J, double convergence_threshold = 1.0e-08, size_t maximum_number_of_iterations = 128);


    // PUBLIC METHODS

    /**
     *  Find a solution to the problem f(x) = 0, where f is a vector function and x is a vector
     */
    void solve(VectorX<double>& x);
};


}  // namespace GQCP
