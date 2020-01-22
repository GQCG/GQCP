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


#include "Mathematical/Algorithm/IterationStep.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationEnvironment.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"

#include <type_traits>


namespace GQCP {


/**
 *  An iteration step that produces updated variables according to a Newton step.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the variables of the system of equations
 *  @tparam _Environment        the type of the calculation environment
 */
template <typename _Scalar, typename _Environment>
class NewtonStepUpdate :
    public IterationStep<_Environment> {

public:
    using Scalar = _Scalar;
    using Environment = _Environment;
    static_assert(std::is_same<Scalar, typename Environment::Scalar>::value, "The scalar type must match that of the environment");
    static_assert(std::is_base_of<NonLinearEquationEnvironment<Scalar>, Environment>::value, "The environment type must derive from NonLinearEquationEnvironment.");


public:

    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Calculate a new iteration of the variables and add them to the environment.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        const auto& x = environment.variables.back();
        const auto& f = environment.f;
        const auto& J = environment.J;

        // Calculate f(x) and J(x), i.e. the values of the vector field and its Jacobian at the given x
        VectorX<Scalar> f_vector = f(x);
        SquareMatrix<Scalar> J_matrix = J(x);
        const VectorX<Scalar> dx = J_matrix.colPivHouseholderQr().solve(-f_vector);  // the actual Newton step

        environment.variables.push_back(x + dx);
    }
};


}  // namespace GQCP
