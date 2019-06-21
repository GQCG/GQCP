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
#include "math/optimization/step.hpp"

#include "math/SquareMatrix.hpp"


namespace GQCP {


/**
 *  @param x        the current point
 *  @param f        a callable vector function
 *  @param J        the corresponding Jacobian function
 *
 *  @return the Newton step
 *      J(x) p = - f
 */
VectorX<double> newtonStep(const VectorX<double>& x, const VectorFunction& f, const MatrixFunction& J) {

    // Calculate f(x) and J(x), i.e. the values of the vector field and its Jacobian at the given x
    VectorX<double> f_vector = f(x);
    SquareMatrix<double> J_matrix = J(x);

    // Return the actual Newton step
    return J_matrix.colPivHouseholderQr().solve(-f_vector);
}


/**
 *  Calculate the electric polarizability from the linear wave function response
 * 
 *  @param F_p          the electric response force (d^2E/dFdp)
 *  @param response     the linear wave function response
 * 
 *  @return the components of the electric polarizability
 */
Matrix<double, 3, 3> calculateElectricPolarizability(const Matrix<double, Dynamic, 3>& F_p, const Matrix<double, Dynamic, 3>& response) {

    return - F_p.transpose() * response;  // no explicit second-order partial perturbation derivative for electrical response
}


}  // namespace GQCP
