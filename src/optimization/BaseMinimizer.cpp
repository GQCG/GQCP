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
#include "optimization/BaseMinimizer.hpp"



namespace GQCP {


/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given initial guess @param x0 and a @param convergence_threshold
 */
BaseMinimizer::BaseMinimizer(const Eigen::VectorXd& x0, double convergence_threshold) :
    x0(x0),
    convergence_threshold(convergence_threshold)
{}


/*
 *  GETTERS
 */
Eigen::VectorXd BaseMinimizer::get_solution() const {
    if (!this->is_solved) {
        throw std::logic_error("The solution hasn't been found and you are trying to get it.");
    } else {
        return this->x;
    }
}


}  // namespace GQCP
