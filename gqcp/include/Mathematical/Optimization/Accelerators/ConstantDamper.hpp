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


#include <deque>



namespace GQCP {


/**
 *  An accelerator that produces its accelerated subject by a linear combination of the two previous subjects.
 */
template <typename _Subject>
class ConstantDamper {
public:
    using Subject = _Subject;


private:
    double alpha;  // the damping factor

    std::deque<Subject> subjects;


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param alpha            the damping factor
     */
    ConstantDamper(const double alpha) :
        alpha (alpha)
    {
        if ((this->alpha > 1.0) || (this->alpha < 0)) {
            throw std::invalid_argument("ConstantDamper::ConstantDamper(const double): The given damping factor must be between 0.0 and 1.0.");
        }
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Add a subject to the subspace of subjects
     */
    Subject accelerate(const Subject& subject) {

        // Add the new subject to the subject subspace and perform an acceleration step if possible
        this->add_subject(subject);

        if (this->canAccelerate()) {
            Subject accelerated_subject = this->subjects.at(1) * this->factor() - this->subjects.at(0) * (1 - this->factor());
            return accelerated_subject;
        }

        else {  // no acceleration is possible
            return subject;
        }
    }


    void add(const Subject& subject) {
        this->subjects.emplace(subject);

        // Remove the oldest subject if there are enough subjects in the subspace
        if (this->subjects.size() >= 2) {
            this->subjects.pop_front();
        }
    }


    /**
     *  @return if this accelerator can produce an accelerated subject
     */
    bool canAccelerate() const {
        if (this->subjects.size() == 2) {
            return true;
        } else {
            return false;
        }
    }
};


}  // namespace GQCP
