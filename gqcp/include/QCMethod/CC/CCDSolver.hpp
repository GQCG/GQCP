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
#include "QCMethod/CC/CCDIntermediatesUpdate.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"
#include "QCMethod/CC/T2DIIS.hpp"
#include "QCMethod/CC/T2ErrorCalculation.hpp"


namespace GQCP {


/**
 *  A factory class that can construct CCD solvers in an easy way.
 * 
 *  @tparam _Scalar             The scalar type that is used to represent the amplitudes.
 */
template <typename _Scalar>
class CCDSolver {
public:
    using Scalar = _Scalar;


public:
    /*
     *  MARK: Factory methods
     */

    /**
     *  Create a plain CCD solver.
     * 
     *  @param threshold                            The threshold that is used in comparing the amplitudes.
     *  @param maximum_number_of_iterations         The maximum number of iterations the algorithm may perform.
     * 
     *  @return A plain CCD solver that uses the norm of the difference of consecutive amplitudes as a convergence criterion.
     */
    static IterativeAlgorithm<CCSDEnvironment<Scalar>> Plain(const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a plain CCD solver.
        StepCollection<CCSDEnvironment<Scalar>> plain_ccd_cycle {};
        plain_ccd_cycle
            .add(CCDIntermediatesUpdate<Scalar>())
            .add(CCDAmplitudesUpdate<Scalar>())
            .add(CCDEnergyCalculation<Scalar>());


        // Create a convergence criterion on the norm of subsequent T2-amplitudes, which is facilitated by the .norm() API of the T2-amplitudes.
        using T2ConvergenceType = ConsecutiveIteratesNormConvergence<T2Amplitudes<Scalar>, CCSDEnvironment<Scalar>>;
        const auto t2_extractor = [](const CCSDEnvironment<Scalar>& environment) { return environment.t2_amplitudes; };
        const T2ConvergenceType t2_convergence_criterion {threshold, t2_extractor, "the T2 amplitudes"};

        // Put together the pieces of the algorithm.
        return IterativeAlgorithm<CCSDEnvironment<Scalar>>(plain_ccd_cycle, t2_convergence_criterion, maximum_number_of_iterations);
    }


    /**
     *  Create a DIIS CCD solver.
     * 
     *  @param minimum_subspace_dimension           The minimum number of T2 amplitudes that have to be in the subspace before enabling DIIS.
     *  @param maximum_subspace_dimension           The maximum number of T2 amplitudes that can be handled by DIIS.
     *  @param threshold                            The threshold that is used in comparing the amplitudes.
     *  @param maximum_number_of_iterations         The maximum number of iterations the algorithm may perform.
     * 
     *  @return A DIIS CCD solver that uses the norm of the difference of consecutive amplitudes as a convergence criterion.
     */
    static IterativeAlgorithm<CCSDEnvironment<Scalar>> DIIS(const size_t minimum_subspace_dimension = 6, const size_t maximum_subspace_dimension = 6, const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a DIIS CCD solver.
        StepCollection<CCSDEnvironment<Scalar>> diis_ccd_cycle {};
        diis_ccd_cycle
            .add(CCDIntermediatesUpdate<Scalar>())
            .add(CCDAmplitudesUpdate<Scalar>())
            .add(T2ErrorCalculation<Scalar>())
            .add(T2DIIS<Scalar>(minimum_subspace_dimension, maximum_subspace_dimension))
            .add(CCDEnergyCalculation<Scalar>());

        // Create a convergence criterion on the norm of subsequent T2-amplitudes, which is facilitated by the .norm() API of the T2-amplitudes.
        using T2ConvergenceType = ConsecutiveIteratesNormConvergence<T2Amplitudes<Scalar>, CCSDEnvironment<Scalar>>;
        const auto t2_extractor = [](const CCSDEnvironment<Scalar>& environment) { return environment.t2_amplitudes; };
        const T2ConvergenceType t2_convergence_criterion {threshold, t2_extractor, "the T2 amplitudes"};

        // Put together the pieces of the algorithm.
        return IterativeAlgorithm<CCSDEnvironment<Scalar>>(diis_ccd_cycle, t2_convergence_criterion, maximum_number_of_iterations);
    }
};


}  // namespace GQCP
