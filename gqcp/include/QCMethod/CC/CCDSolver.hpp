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


#include "Mathematical/Algorithm/ConvergenceCriterion.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Mathematical/Algorithm/StepCollection.hpp"
#include "Mathematical/Optimization/ConsecutiveIteratesNormConvergence.hpp"
#include "QCMethod/CC/CCDAmplitudesUpdate.hpp"
#include "QCMethod/CC/CCDEnergyCalculation.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"
#include "QCMethod/CC/CCDIntermediatesUpdate.hpp"


namespace GQCP {


/**
 *  A factory class that can construct CCD solvers in an easy way.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the amplitudes
 */
template <typename _Scalar>
class CCDSolver {
public:
    using Scalar = _Scalar;


public:
    /*
     *  PUBLIC STATIC METHODS
     */

    /**
     *  @param threshold                            the threshold that is used in comparing the amplitudes
     *  @param maximum_number_of_iterations         the maximum number of iterations the algorithm may perform
     * 
     *  @return a plain CCD solver that uses the norm of the difference of consecutive amplitudes as a convergence criterion
     */
    static IterativeAlgorithm<CCSDEnvironment<Scalar>> Plain(const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a plain CCD solver.
        StepCollection<CCSDEnvironment<Scalar>> plain_ccd_cycle {};
        plain_ccd_cycle
            .add(CCDIntermediatesUpdate<Scalar>())
            .add(CCDAmplitudesUpdate<Scalar>())
            .add(CCDEnergyCalculation<Scalar>());


        // Create a compound convergence criterion on the norm of subsequent T2-amplitudes, which is facilitated by the .norm() API of the T2-amplitudes.
        using T2ConvergenceType = ConsecutiveIteratesNormConvergence<T2Amplitudes<Scalar>, CCSDEnvironment<Scalar>>;
        const auto t2_extractor = [](const CCSDEnvironment<Scalar>& environment) { return environment.t2_amplitudes; };
        const T2ConvergenceType t2_convergence_criterion {threshold, t2_extractor, "the T2 amplitudes"};


        // Put together the pieces of the algorithm.
        return IterativeAlgorithm<CCSDEnvironment<Scalar>>(plain_ccd_cycle, t2_convergence_criterion, maximum_number_of_iterations);
    }
};


}  // namespace GQCP
