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
#include "step.hpp"


namespace numopt {


/**
 *  @param x        the current point
 *  @param f        a callable vector function
 *  @param J        the corresponding Jacobian function
 *
 *  @return the Newton step
 *      J(x) p = - f
 */
Eigen::VectorXd newtonStep(const Eigen::VectorXd& x, const VectorFunction& f, const MatrixFunction& J) {

    // Calculate f(x) and J(x), i.e. the values of the vector field and its Jacobian at the given x
    Eigen::VectorXd f_vector = f(x);
    Eigen::MatrixXd J_matrix = J(x);

    // Return the actual Newton step
    return J_matrix.colPivHouseholderQr().solve(-f_vector);
}


}  // namespace numopt
