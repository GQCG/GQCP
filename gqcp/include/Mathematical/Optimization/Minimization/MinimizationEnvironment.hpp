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
 *  An environment that can be used to minimize scalar functions.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the variables of scalar function
 */
template <typename _Scalar>
class MinimizationEnvironment:
    public OptimizationEnvironment<VectorX<_Scalar>> {

public:
    using Scalar = _Scalar;


public:
    VectorFunction<Scalar> gradient_function;  // a callable function that produces the gradient of the scalar function, evaluated at the given variables
    MatrixFunction<Scalar> hessian_function;   // a callable function that produces the Hessian of the scalar function, evaluated at the given variables

    std::deque<double> function_values;  // values for the evaluated scalar function (often a sort of 'cost' function)


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize the optimization environment with an initial guess
     * 
     *  @param initial_guess                the initial guess for the variables
     *  @param gradient_function            a callable function that produces the gradient of the scalar function, evaluated at the given variables
     *  @param hessian_function             a callable function that produces the Hessian of the scalar function, evaluated at the given variables
     */
    MinimizationEnvironment(const VectorX<_Scalar>& initial_guess, const VectorFunction<Scalar>& gradient_function, const MatrixFunction<Scalar>& hessian_function) :
        OptimizationEnvironment<VectorX<_Scalar>>(initial_guess),
        gradient_function {gradient_function},
        hessian_function {hessian_function} {}
};


}  // namespace GQCP
