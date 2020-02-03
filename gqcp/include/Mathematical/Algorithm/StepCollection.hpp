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
#include "Utilities/memory.hpp"
#include "Utilities/type_traits.hpp"

#include <vector>


namespace GQCP {


/**
 *  A collection of steps to be executed in a consecutive order.
 * 
 *  This iteration cycle maintains the ownership of its constituting steps.
 * 
 *  @param _Environment             the type of the environment that this iteration step can read from and write to
 */
template <typename _Environment>
class StepCollection :
    public Step<_Environment> {
public:
    using Environment = _Environment;


private:
    std::vector<std::shared_ptr<Step<Environment>>> steps;  // the consecutive steps that this collection consists of


public:

    /*
     *  PUBLIC METHODS
     */

    /**
     *  Add a new step to the collection of steps.
     * 
     *  @return the modified collection of steps, in order to allow chaining.
     */
    template <typename Z = Step<Environment>>
    enable_if_t<std::is_same<Environment, typename Z::Environment>::value, StepCollection<Environment>&> add(const Z& step) {
        this->steps.push_back(std::make_shared<Z>(step));
        return *this;
    }


    /**
     *  Execute all the steps in this collection.
     * 
     *  @param environment              the environment that this step can read from and write to
     */
    void execute(Environment& environment) {
        for (const auto& step : this->steps) {
            step->execute(environment);
        }
    }
};


}  // namespace GQCP
