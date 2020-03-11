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


#include "Mathematical/Optimization/OptimizationEnvironment.hpp"
#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/**
 *  An environment that can be used to solve linear systems of equations.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the variables of the system of equations
 */
template <typename _Scalar>
class LinearEquationEnvironment {
public:
    using Scalar = _Scalar;


public:
    MatrixX<Scalar> A;  // the matrix that corresponds to the left-hand side of the linear system of equations
    MatrixX<Scalar> b;  // the matrix/vector that corresponds to the right-hand side of the linear system of equations

    MatrixX<Scalar> x;  // the matrix/vector of solutions


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param A            the matrix that corresponds to the left-hand side of the linear system of equations
     *  @param b            the matrix/vector that corresponds to the right-hand side of the linear system of equations
     */
    LinearEquationEnvironment(const MatrixX<Scalar>& A, const MatrixX<Scalar>& b) :
        A (A),
        b (b)
    {}
};


}  // namespace GQCP
