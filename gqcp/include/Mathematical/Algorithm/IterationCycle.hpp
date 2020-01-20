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
#include "Utilities/memory.hpp"
#include "Utilities/type_traits.hpp"

#include <vector>


namespace GQCP {


/**
 *  A collection of iteration steps that constitutes one iteration. An iteration is defined as the part of an iterative algorithm that happens between convergence checks.
 * 
 *  This iteration cycle maintains the ownership of its constituting steps.
 * 
 *  @param _Environment             the type of the environment that this iteration step can read from and write to
 */
template <typename _Environment>
class IterationCycle {
public:
    using Environment = _Environment;


private:
    std::vector<std::unique_ptr<IterationStep<Environment>>> steps;  // the consecutive steps that this iteration cycle consists of


public:

    /*
     *  CONSTRUCTORS
     */
    IterationCycle() : {}


    /**
     *  Add a new step to the iteration cycle.
     * 
     *  @return the modified iteration cycle, in order to allow chaining.
     */
    template <typename Z = IterationStep>
    std::enable_if_t<std::is_same<Environment, typename Z::Environment>, IterationCycle&> add(const Z<Environment>& step) {
        this->steps.append(std::make_unique(step));
        return *this;
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Execute all the steps of this iteration cycle.
     * 
     *  @param environment              the environment that this iteration step can read from and write to
     */
    void execute(Environment& environment) {
        for (const auto& step : this->steps) {
            step->execute(environment);
        }
    }
};


}  // namespace GQCP
