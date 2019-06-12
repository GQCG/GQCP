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
#ifndef GQCP_ORBITALOPTIMIZATIONOPTIONS_HPP
#define GQCP_ORBITALOPTIMIZATIONOPTIONS_HPP


#include <cstdlib>


namespace GQCP {


/**
 *  A class that holds options for orbital optimization
 */
class OrbitalOptimizationOptions {
private:
    bool should_minimize;  // if the algorithm should look for a minimum or not
    double convergence_threshold;
    size_t maximum_number_of_iterations;


public:
    // CONSTRUCTORS

    /**
     *  @param should_minimize                  if the algorithm should look for a minimum or not
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     */
    OrbitalOptimizationOptions(const bool should_minimize=true, const double convergence_threshold = 1.0e-08, const double maximum_number_of_iterations = 128);


    // NAMED CONSTRUCTORS

    /**
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     * 
     *  @return orbital optimization options suited for the maximization of a cost function
     */
    static OrbitalOptimizationOptions OrbitalMaximizationOptions(const double convergence_threshold = 1.0e-08, const double maximum_number_of_iterations = 128);

    /**
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     * 
     *  @return orbital optimization options suited for the minimization of a cost function
     */
    static OrbitalOptimizationOptions OrbitalMinimizationOptions(const double convergence_threshold = 1.0e-08, const double maximum_number_of_iterations = 128);


    // PUBLIC METHODS

    /**
     *  @return if the algorithm should look for a minimum or not
     */
    bool shouldMinimize() const { return this->should_minimize; }

    /**
     *  @return if the algorithm should look for a maximum or not
     */
    bool shouldMaximize() const { return !this->shouldMinimize(); }

    /**
     *  @return the threshold used to check for convergence
     */
    double convergenceThreshold() const { return this->convergence_threshold; }

    /**
     *  @return the maximum number of iterations that may be used to achieve convergence
     */
    size_t maximumNumberOfIterations() const { return this->maximum_number_of_iterations; }
};


}  // namespace GQCP


#endif  // GQCP_ORBITALOPTIMIZATIONOPTIONS_HPP
