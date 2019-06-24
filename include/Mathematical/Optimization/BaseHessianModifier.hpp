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
#ifndef GQCP_BASEHESSIANMODIFIER_HPP
#define GQCP_BASEHESSIANMODIFIER_HPP


#include "Mathematical/SquareMatrix.hpp"


namespace GQCP {


/**
 *  A base functor for Hessian modifiers
 */
class BaseHessianModifier {
public:
    // DESTRUCTOR
    virtual ~BaseHessianModifier() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  @param hessian      the current indefinite Hessian
     * 
     *  @return a modified Hessian that is made positive (for minimizers) or negative (for maximizers) definite
     */
    virtual SquareMatrix<double> operator()(const SquareMatrix<double>& hessian) = 0;
};

/**
 *  Calculate the electric polarizability from the linear wave function response
 * 
 *  @param F_p          the electric response force (d^2E/dFdp)
 *  @param response     the linear wave function response
 * 
 *  @return the components of the electric polarizability
 */
Matrix<double, 3, 3> calculateElectricPolarizability(const Matrix<double, Dynamic, 3>& F_p, const Matrix<double, Dynamic, 3>& response);


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
Matrix<double, 3, 3> calculateElectricPolarizability(const Matrix<double, Dynamic, 3>& F_p, const Matrix<double, Dynamic, 3>& x, const Matrix<double, Dynamic, 3>& A_lambda, const Matrix<double, Dynamic, 3>& y);


}  // namespace GQCP



#endif  // GQCP_BASEHESSIANMODIFIER_HPP
