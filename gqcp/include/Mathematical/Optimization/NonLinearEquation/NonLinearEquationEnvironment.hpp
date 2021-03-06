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

#pragma once


#include "Mathematical/Optimization/OptimizationEnvironment.hpp"
#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/**
 *  An environment that can be used to solve non-linear systems of equations
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the variables of the system of equations
 */
template <typename _Scalar>
class NonLinearEquationEnvironment:
    public OptimizationEnvironment<VectorX<_Scalar>> {

public:
    using Scalar = _Scalar;


public:
    VectorFunction<Scalar> f;  // a callable function that produces a vector function that represents the system of equations at the given variables
    MatrixFunction<Scalar> J;  // a callable function that produces a matrix that represents the Jacobian of the system of equations at the given variables


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize the optimization environment with an initial guess
     * 
     *  @param initial_guess                the initial guess for the variables
     *  @param f                            a callable function that produces a vector function that represents the system of equations at the given variables
     *  @param J                            a callable function that produces a matrix that represents the Jacobian of the system of equations at the given variables
     */
    NonLinearEquationEnvironment(const VectorX<_Scalar>& initial_guess, const VectorFunction<Scalar>& f, const MatrixFunction<Scalar>& J) :
        OptimizationEnvironment<VectorX<_Scalar>>(initial_guess),
        f {f},
        J {J} {}
};


}  // namespace GQCP
