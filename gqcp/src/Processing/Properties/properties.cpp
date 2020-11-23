// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "Processing/Properties/properties.hpp"

#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "Processing/Properties/expectation_values.hpp"


namespace GQCP {


/**
 *  Calculate the electric polarizability from the linear wave function response
 * 
 *  @param F_p          the electric response force (d^2E/dFdp)
 *  @param x            the linear wave function response
 * 
 *  @return the components of the electric polarizability
 */
Matrix<double, 3, 3> calculateElectricPolarizability(const Matrix<double, Dynamic, 3>& F_p, const Matrix<double, Dynamic, 3>& x) {

    // No explicit second-order partial perturbation derivative for electrical response
    return -(x.transpose() * F_p);  // minus sign because of definition of electric polarizability
}


/**
 *  Calculate the electric polarizability from the linear wave function response and the linear multiplier response
 * 
 *  @param F_p              the electric parameter response force
 *  @param x                the linear wave function response
 *  @param A_lambda         the first part of the electric multiplier response force
 *  @param y                the linear multiplier response
 * 
 *  @return the components of the electric polarizability
 */
Matrix<double, 3, 3> calculateElectricPolarizability(const Matrix<double, Dynamic, 3>& F_p, const Matrix<double, Dynamic, 3>& x, const Matrix<double, Dynamic, 3>& A_lambda, const Matrix<double, Dynamic, 3>& y) {

    // No explicit second-order partial perturbation derivative for electrical response
    return -(x.transpose() * A_lambda + y.transpose() * F_p);  // minus sign because of definition of electric polarizability
}

}  // namespace GQCP
