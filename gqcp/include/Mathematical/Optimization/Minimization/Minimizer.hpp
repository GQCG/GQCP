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


#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Mathematical/Optimization/ConsecutiveIteratesNormConvergence.hpp"
#include "Mathematical/Optimization/Minimization/MinimizationEnvironment.hpp"
#include "Mathematical/Optimization/Minimization/NewtonStepUpdate.hpp"
#include "Mathematical/Optimization/OptimizationEnvironment.hpp"


namespace GQCP {


/**
 *  A factory class that can construct minimizers in an easy way.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the variables of function that is minimized
 */
template <typename _Scalar>
class Minimizer {
public:
    using Scalar = _Scalar;


public:
    /*
     * STATIC PUBLIC METHODS
     */

    /**
     *  @param threshold                            the threshold that is used in comparing the density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the algorithm may perform
     * 
     *  @return a Newton-step minimizer that uses the norm of the difference of two consecutive iterations of variables as a convergence criterion
     */
    static IterativeAlgorithm<MinimizationEnvironment<Scalar>> Newton(const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a Newton-step minimizer
        StepCollection<MinimizationEnvironment<Scalar>> newton_cycle {};
        newton_cycle.add(GQCP::Minimization::NewtonStepUpdate<Scalar, MinimizationEnvironment<Scalar>>());

        // Create a convergence criterion on the norm of subsequent iterations of variables
        const ConsecutiveIteratesNormConvergence<VectorX<Scalar>, MinimizationEnvironment<Scalar>> convergence_criterion {threshold};

        return IterativeAlgorithm<MinimizationEnvironment<Scalar>>(newton_cycle, convergence_criterion, maximum_number_of_iterations);
    }
};


}  // namespace GQCP
