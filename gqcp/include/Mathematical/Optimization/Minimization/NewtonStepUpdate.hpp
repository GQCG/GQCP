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


#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/Minimization/MinimizationEnvironment.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NewtonStepUpdate.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationEnvironment.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"

#include <type_traits>


namespace GQCP {
namespace Minimization {


/**
 *  An iteration step that produces updated variables according to a Newton step.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the variables of the system of equations
 *  @tparam _Environment        the type of the calculation environment
 */
template <typename _Scalar, typename _Environment>
class NewtonStepUpdate:
    public Step<_Environment> {

public:
    using Scalar = _Scalar;
    using Environment = _Environment;
    static_assert(std::is_same<Scalar, typename Environment::Scalar>::value, "The scalar type must match that of the environment");
    static_assert(std::is_base_of<MinimizationEnvironment<Scalar>, Environment>::value, "The environment type must derive from MinimizationEnvironment.");


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

        // Use the non-linear equation Newton step as an implementation
        const auto& variables = environment.variables.back();
        const auto& f = environment.gradient_function;
        const auto& J = environment.hessian_function;

        GQCP::NonLinearEquationEnvironment<Scalar> non_linear_environment {variables, f, J};
        GQCP::NonLinearEquation::NewtonStepUpdate<Scalar, NonLinearEquationEnvironment<Scalar>>().execute(non_linear_environment);  // this adds the Newton-step updated variables to the non-linear environment

        const auto& new_variables = non_linear_environment.variables.back();
        environment.variables.push_back(new_variables);
    }
};


}  // namespace Minimization
}  // namespace GQCP
