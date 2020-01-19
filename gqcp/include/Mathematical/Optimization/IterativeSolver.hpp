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


#include <cstddef>


namespace GQCP {


/**
 *  A base class for iterative solvers: at every step, the iterate is updated.
 * 
 *  @tparam _Iterate            the type of the iterate
 */
template <typename _Iterate>
class IterativeSolver {
public:
    using Iterate = _Iterate;


private:
    size_t maximum_number_of_iterations;
    size_t iteration = 0;  // the current iteration counter


protected:
    Iterate iterate;  // the iterate, which is continuously updated in every iteration step


public:

    /*
     *  CONSTRUCTORS
     */
    
    /**
     *  Initialize the solver with an initial guess
     * 
     *  @param initial_guess                        the initial guess to the solver
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform
     */
    IterativeSolver(const Iterate& initial_guess, const size_t maximum_number_of_iterations = 128) :
        maximum_number_of_iterations (maximum_number_of_iterations),
        iterate (initial_guess)
    {}


    /*
     *  PUBLIC PURE VIRTUAL METHODS
     */

    /**
     *  @return if the algorithm is considered to be converged
     */
    virtual bool isConverged() = 0;

    /**
     *  @return a new iterate to be used in the next iteration
     */
    virtual Iterate calculateNextIterate() = 0;



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the maximum number of iterations this solver is allowed to perform
     */
    size_t maximumNumberOfIterations() const { return this->maximum_number_of_iterations; }

    /**
     *  @return the number of iterations that this solver has performed
     */
    size_t numberOfIterations() const { return this->iteration; }

    /**
     *  Iterate until the algorithm has converged, i.e. the actual solving step
     * 
     *  @return the solution of the algorithm
     */
    Iterate solve() {

        while (!this->isConverged()) {
            this->iterate = this->calculateNextIterate();

            this->iteration++;
            if (this->iteration >= this->maximum_number_of_iterations) {
                throw std::runtime_error("IterativeSolver::solve(): The iterative solver did not converge in the given maximum number of iterations.");
            }
        }

        // If the algorithm has converged, we are safe to return the final iterate
        return this->iterate;
    }
};


}  // namespace GQCP
