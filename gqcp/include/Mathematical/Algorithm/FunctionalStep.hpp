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


#include "Mathematical/Algorithm/Step.hpp"

#include <functional>

namespace GQCP {


/**
 *  A general algorithmic step that acts as a wrapper around a function call.
 * 
 *  @param _Environment             the type of the environment that this step can read from and write to
 */
template <typename _Environment>
class FunctionalStep:
    public Step<_Environment> {

public:
    using Environment = _Environment;
    using Function = std::function<void(Environment&)>;


private:
    Function function;  // the function that this Step wraps


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Create a Step by wrapping a function.
     * 
     *  @param function             the function that this Step wraps
     */
    FunctionalStep(const Function& function) :
        function {function} {}


    /*
     *  PUBLIC OVERRIDDEN FUNCTIONS
     */

    /**
     *  Execute/perform this algorithm step.
     * 
     *  @param environment              the environment that this step can read from and write to
     */
    void execute(Environment& environment) override {
        this->function(environment);
    }
};


}  // namespace GQCP
