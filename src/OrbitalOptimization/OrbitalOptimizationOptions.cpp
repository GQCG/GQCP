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

#include "OrbitalOptimization/OrbitalOptimizationOptions.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param should_minimize                  if the algorithm should look for a minimum or not
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
 */
OrbitalOptimizationOptions::OrbitalOptimizationOptions(const bool should_minimize, const double convergence_threshold, const double maximum_number_of_iterations) :
    should_minimize (should_minimize),
    convergence_threshold (convergence_threshold),
    maximum_number_of_iterations (maximum_number_of_iterations)
{}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
 * 
 *  @return orbital optimization options suited for the maximization of a cost function
 */
OrbitalOptimizationOptions OrbitalOptimizationOptions::OrbitalMaximizationOptions(const double convergence_threshold, const double maximum_number_of_iterations) {
    return OrbitalOptimizationOptions(false, convergence_threshold, maximum_number_of_iterations);
}

/**
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
 * 
 *  @return orbital optimization options suited for the minimization of a cost function
 */
OrbitalOptimizationOptions OrbitalOptimizationOptions::OrbitalMinimizationOptions(const double convergence_threshold, const double maximum_number_of_iterations) {
    return OrbitalOptimizationOptions(true, convergence_threshold, maximum_number_of_iterations);
}


}  // namespace GQCP