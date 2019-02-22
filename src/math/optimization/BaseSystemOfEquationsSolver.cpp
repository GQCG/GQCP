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
#include "math/optimization/BaseSystemOfEquationsSolver.hpp"



namespace GQCP {


/*
 *  CONSTRUCTORS
 */
/**
 *  @param x0                               an initial guess
 *  @param convergence_threshold            the threshold for convergence on the norm of the gradient
 *  @param maximum_number_of_iterations     the maximum number of iterations in the algorithm
 */
BaseSystemOfEquationsSolver::BaseSystemOfEquationsSolver(const Eigen::VectorXd& x0, double convergence_threshold, size_t maximum_number_of_iterations) :
    x (x0),
    convergence_threshold (convergence_threshold),
    maximum_number_of_iterations (maximum_number_of_iterations)
{}



/*
 *  GETTERS
 */
const Eigen::VectorXd& BaseSystemOfEquationsSolver::get_solution() const {
    if (!this->is_solved) {
        throw std::logic_error("The solution hasn't been found and you are trying to get it.");
    } else {
        return this->x;
    }
}


}  // namespace GQCP
