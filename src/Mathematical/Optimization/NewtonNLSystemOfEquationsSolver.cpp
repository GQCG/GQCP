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
#include "Mathematical/Optimization/NewtonNLSystemOfEquationsSolver.hpp"

#include "Mathematical/Optimization/step.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param f                                a callable vector function that returns the value of each of the equations in a resulting vector, given the current variables
 *  @param J                                the corresponding callable Jacobian
 *  @param convergence_threshold            the threshold on the norm of the update that determines convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations in the algorithm
 */
NewtonNLSystemOfEquationsSolver::NewtonNLSystemOfEquationsSolver(const VectorFunction& f, const MatrixFunction& J, double convergence_threshold, size_t maximum_number_of_iterations) :
    convergence_threshold (convergence_threshold), 
    maximum_number_of_iterations (maximum_number_of_iterations),
    f (f),
    J (J)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Find a solution to the problem f(x) = 0, where f is a vector function and x is a vector
 */
void NewtonNLSystemOfEquationsSolver::solve(VectorX<double>& x) {

    size_t iteration_counter = 0;
    bool is_solved = false;
    while (!is_solved) {


        std::cout << "x: " << std::endl << x << std::endl << std::endl;


        // Calculate the Newton step and update the current variables
        const VectorX<double> dx = newtonStep(x, this->f, this->J);
        x += dx;


        // Check for convergence on the norm of the difference
        if (dx.norm() <= this->convergence_threshold) {
            is_solved = true;
        } else {  // not solved yet
            iteration_counter++;

            // If we reach more than the maximum number of iterations, the system is considered not to be converging
            if (iteration_counter >= this->maximum_number_of_iterations) {
                throw std::runtime_error("NewtonNLSystemOfEquationsSolver::solve(VectorX<double>&): The Newton procedure did not converge.");
            }
        }
    }  // while loop
}


}  // namespace GQCP
