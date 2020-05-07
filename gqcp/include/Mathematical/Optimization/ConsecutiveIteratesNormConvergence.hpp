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

#include <deque>
#include <functional>
#include <type_traits>


namespace GQCP {


/**
 *  A convergence criterion that checks if the norm of the difference of two iterates is converged.
 * 
 *  @tparam _Iterate            the type of the iterative variables
 *  @tparam _Environment        the type of the calculation environment
 */
template <typename _Iterate, typename _Environment>
class ConsecutiveIteratesNormConvergence:
    public ConvergenceCriterion<_Environment> {

public:
    using Iterate = _Iterate;
    using Scalar = typename Iterate::Scalar;
    using Environment = _Environment;
    static_assert(std::is_same<Scalar, typename Environment::Scalar>::value, "The scalar types of the iterate and environment must match.");


private:
    double threshold;  // the threshold that is used in comparing the iterates

    std::string iterate_description;  // the description of the the iterates that are compared

    std::function<std::deque<Iterate>(const Environment&)> extractor;  // a function that can extract the correct iterates from the environment, as it's not mandatory to check convergence on the variables, but any iterate (whose .norm() can be calculated) can in principle be used


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param threshold                    the threshold that is used in comparing the iterates
     *  @param extractor                    a function that can extract the correct iterates from the environment. The default is to check the environment on a property called 'variables'
     *  @param iterate_description          the description of the the iterates that are compared
     */
    ConsecutiveIteratesNormConvergence(
        const double threshold = 1.0e-08, const std::function<std::deque<Iterate>(const Environment&)> extractor = [](const Environment& environment) { return environment.variables; }, const std::string& iterate_description = "a general iterate") :
        threshold {threshold},
        extractor {extractor},
        iterate_description {iterate_description} {}


    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return (boost::format("A convergence criterion that checks if the norm of the difference of two iterates (%1%) is converged") % this->iterate_description).str();
    }


    /**
     *  @param environment                  the environment that acts as a sort of calculation space
     * 
     *  @return if the difference of the two most recent iterates has a zero norm, within the tolerance
     */
    bool isFulfilled(Environment& environment) override {

        const auto iterates = this->extractor(environment);

        if (iterates.size() < 2) {
            return false;  // we can't calculate convergence
        }

        // Get the two most recent density matrices and compare the norm of their difference
        const auto second_to_last_it = iterates.end() - 2;  // 'it' for 'iterator'
        const auto& previous = *second_to_last_it;          // dereference the iterator
        const auto current = iterates.back();

        return ((current - previous).norm() <= this->threshold);
    }
};


}  // namespace GQCP
