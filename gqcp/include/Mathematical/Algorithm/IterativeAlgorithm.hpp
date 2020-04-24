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
#include "Mathematical/Algorithm/StepCollection.hpp"
#include <cstddef>


namespace GQCP {


/**
 *  An algorithm that performs iterations. In every iteration, convergence is checked and a set of iteration steps is performed.
 * 
 *  @param _Environment             the type of environment that this algorithm is associated to
 */
template <typename _Environment>
class IterativeAlgorithm {
public:
    using Environment = _Environment;


private:
    size_t maximum_number_of_iterations;
    size_t iteration = 0;  // the number of iterations that have been performed

    StepCollection<Environment> steps;  // the collection of algorithm steps that is performed in-between convergence checks
    std::shared_ptr<ConvergenceCriterion<Environment>> convergence_criterion;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize the members of the iterative algorithm
     * 
     *  @tparam Criterion                           the type of the convergence criterion that is used
     * 
     *  @param steps                                the collection of algorithm steps that is performed in-between convergence checks
     *  @param convergence_criterion                the convergence criterion that must be fulfilled in order for the algorithm to have converged
     *  @param maximum_number_of_iterations         the maximum number of iterations the algorithm may perform
     */
    template <typename Criterion>
    IterativeAlgorithm(const StepCollection<Environment>& steps, const Criterion& convergence_criterion, const size_t maximum_number_of_iterations = 128) :
        maximum_number_of_iterations {maximum_number_of_iterations},
        steps {steps},
        convergence_criterion {std::make_shared<Criterion>(convergence_criterion)} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Perform the iteration steps until convergence is achieved
     * 
     *  @param environment                          the environment that this algorithm is associated to
     */
    void perform(Environment& environment) {

        for (this->iteration = 0; this->iteration <= this->maximum_number_of_iterations; this->iteration++) {  // do at maximum the maximum allowed number of iterations

            // Every iteration consists of two parts:
            //      - the convergence check, which checks if the iterations may stop
            //      - the iteration cycle, i.e. what happens in-between the convergence checks
            if (this->convergence_criterion->isFulfilled(environment)) {
                return;  // exit the loop and function early
            }

            this->steps.execute(environment);
        }

        // Since we will exit the function early if convergence is achieved, the algorithm is considered non-converging if the loop is done.
        throw std::runtime_error("IterativeAlgorithm<Environment>::perform(Environment&): The algorithm didn't find a solution within the maximum number of iterations.");
    }


    /**
     *  @return the number of iterations that have been performed
     */
    size_t numberOfIterations() const { return this->iteration; }
};


}  // namespace GQCP
