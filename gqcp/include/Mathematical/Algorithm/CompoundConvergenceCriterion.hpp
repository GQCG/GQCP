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

#include <boost/format.hpp>

#include <memory>
#include <vector>


namespace GQCP {


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
    using Self = CompoundConvergenceCriterion<Environment>;


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


    /**
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {

        std::string description_string = (boost::format("A compound convergence criterion step consisting of %1% combined criteria:\n") % this->count()).str();

        for (size_t i = 0; i < this->count(); i++) {
            const auto& criterion = this->criteria[i];

            description_string += "\t";
            description_string += std::to_string(i + 1);  // +1 because of computer indices
            description_string += ". ";
            description_string += criterion->description();
            description_string += "\n";
        }
        return description_string;
    }


    /**
     *  @param environment              the environment that this criterion can read from
     * 
     *  @return if this criterion is fulfilled
     */
    bool isFulfilled(Environment& environment) override {

        // Check if every criterion is fulfilled.
        for (const auto& criterion : this->criteria) {
            if (!criterion->isFulfilled(environment)) {
                return false;
            }
        }
        return true;
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Add another convergence criterion that should be fulfilled in this compound convergence criterion.
     * 
     *  @param criterion                    the new convergence criterion
     *  
     *  @return an updated version of *this, in order to allow chaining
     */
    template <typename C>
    Self add(const C& criterion) {
        this->criteria.push_back(std::make_shared<C>(criterion));
        return *this;
    }


    /**
     *  @return the number of simple convergence criteria that are compounded in this one
     */
    size_t count() const { return this->criteria.size(); }
};


}  // namespace GQCP
