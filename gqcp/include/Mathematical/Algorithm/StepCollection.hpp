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
#include "Utilities/memory.hpp"
#include "Utilities/type_traits.hpp"

#include <boost/format.hpp>

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
class StepCollection: public Step<_Environment> {
public:
    using Environment = _Environment;


private:
    std::vector<std::shared_ptr<Step<Environment>>> steps;  // the consecutive steps that this collection consists of


public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {

        std::string description_string = (boost::format("An algorithmic step consisting of %s algorithmic steps:\n") % this->numberOfSteps()).str();

        for (size_t i = 0; i < this->numberOfSteps(); i++) {
            const auto& step = this->steps[i];
            description_string += (boost::format("\t%s. %s\n") % std::to_string(i + 1) % step->description()).str();
        }
        return description_string;
    }


    /**
     *  Execute all the steps in this collection.
     * 
     *  @param environment              the environment that this step can read from and write to
     */
    void execute(Environment& environment) override {
        for (const auto& step : this->steps) {
            step->execute(environment);
        }
    }


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
     *  Insert an algorithm step at the given index.
     * 
     *  @param step                 the step that should be inserted into this algorithm step collection
     *  @param index                the zero-based index that the given step should be performed at in the total collection of steps
     */
    template <typename Z = Step<Environment>>
    enable_if_t<std::is_same<Environment, typename Z::Environment>::value, void> insert(const Z& step, const size_t index) {

        // Check if the index is out of bounds.
        if (index > this->numberOfSteps()) {
            throw std::invalid_argument("StepCollection::insert(const Z&, const size_t): Cannot insert at the given index.");
        }


        // Actually insert in the std::vector of steps.
        if (index < this->numberOfSteps()) {  // within the current bounds
            const auto it = this->steps.begin();
            this->steps.insert(it + index, std::make_shared<Z>(step));  // inserting at an given index goes through an iterator
        }

        else if (index == this->numberOfSteps()) {  // at the end
            this->add(step);
        }
    }


    /**
     *  @return the number of steps that are in this consecutive collection
     */
    size_t numberOfSteps() const { return this->steps.size(); }


    /**
     *  Remove the step at the given index.
     * 
     *  @param index                the zero-based index of the step in this collection that should be removed
     */
    void remove(const size_t index) {
        this->steps.erase(this->steps.begin() + index);
    }


    /**
     *  Replace an algorithm step at the given index.
     * 
     *  @param step                 the step that should be inserted into this algorithm
     *  @param index                the zero-based index of the step that should be replaced
     */
    template <typename Z = Step<Environment>>
    enable_if_t<std::is_same<Environment, typename Z::Environment>::value, void> replace(const Z& step, const size_t index) {
        this->steps[index] = std::make_shared<Z>(step);
    }
};


}  // namespace GQCP
