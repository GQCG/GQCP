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


#include "Mathematical/Algorithm/StepCollection.hpp"


namespace GQCP {


/**
 *  An algorithm that only performs one collection of steps.
 * 
 *  @param _Environment             the type of environment that this algorithm is associated to
 */
template <typename _Environment>
class Algorithm {
public:
    using Environment = _Environment;


private:
    StepCollection<Environment> steps;  // the steps that should be executed when this algorithm is performed


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param steps                      the steps that should be executed when this algorithm is performed
     */
    Algorithm(const StepCollection<Environment>& steps) :
        steps {steps} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return a textual description of this iterative algorithm.
     */
    std::string description() const {

        std::string description_string = "An algorithm consisting of the following steps:\n";
        description_string += this->steps.description();
        return description_string;
    }


    /**
     *  @param environment              the environment that this algorithm is associated to
     */
    void perform(Environment& environment) {
        this->steps.execute(environment);
    }
};


}  // namespace GQCP
