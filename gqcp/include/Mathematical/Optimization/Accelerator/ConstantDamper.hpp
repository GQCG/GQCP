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


namespace GQCP {


/**
 *  An accelerator that produces its accelerated subject by a damped linear combination of the two previous subjects.
 */
class ConstantDamper {
private:
    double alpha;  // the damping factor


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
     *  Calculate an accelerated subject using a constant damping step.
     * 
     *  @tparam Subject                 the type of subject that should be accelerated
     * 
     *  @param last                     the most recent subject
     *  @param next_to_last             the second most recent subject
     * 
     *  @return an accelerated subject
     */
    template <typename Subject>
    Subject accelerate(const Subject& last, const Subject& next_to_last) {

        Subject accelerated_subject = last * this->alpha + next_to_last * (1 - this->alpha);
        return accelerated_subject;
    }
};


}  // namespace GQCP
