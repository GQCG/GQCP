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

#include <memory>
#include <vector>


namespace GQCP {


/**
 *  FIXME: rewrite this using std::shared_ptr<ConvergenceCriterion>, see IterativeAlgorithm
 */

/**
 *  A criterion that checks convergence by requiring that multiple convergence criteria are fulfilled simultaneously.
 * 
 *  @param _Environment             the type of the environment that this criterion can read from
 */
template <typename _Environment>
class CompoundConvergenceCriterion:
    public ConvergenceCriterion<_Environment> {

public:
    using Environment = _Environment;


private:
    std::vector<std::shared_ptr<ConvergenceCriterion<Environment>>> criteria;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  A basic constructor from two criteria that should be fulfilled at the same time.
     * 
     *  @param criterion1           the first convergence criterion
     *  @param criterion2           the second convergence criterion
     */
    template <typename C1, typename C2>
    CompoundConvergenceCriterion(const C1& criterion1, const C2& criterion2) :
        criteria {std::make_shared<C1>(criterion1), std::make_shared<C2>(criterion2)} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param environment              the environment that this criterion can read from
     * 
     *  @return if this criterion is fulfilled
     */
    bool isFulfilled(Environment& environment) override {
        return this->criterion1.isFulfilled(environment) && this->criterion2.isFulfilled(environment);
    }
};


}  // namespace GQCP
