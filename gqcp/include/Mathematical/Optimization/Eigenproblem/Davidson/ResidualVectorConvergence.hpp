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


#include "Mathematical/Algorithm/ConvergenceCriterion.hpp"



namespace GQCP {


/**
 *  A convergence criterion that checks if the norm of each of the residual vectors is smaller than a threshold.
 * 
 *  @tparam _Environment        the type of the calculation environment
 */
template <typename _Environment>
class ResidualVectorConvergence :
    public ConvergenceCriterion<_Environment> {

public:
    using Environment = _Environment;


private:
    double threshold;  // the threshold that is used in checking the norm of the residual vectors


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param threshold                    the threshold that is used in checking the norm of the residual vectors
     */
    ResidualVectorConvergence(const double threshold = 1.0e-08) :
        threshold (threshold)
    {}


    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @param environment                  the environment that acts as a sort of calculation space
     * 
     *  @return if the norm of each of the residual vectors is smaller than a threshold
     */
    bool isFulfilled(Environment& environment) override {

        const auto R = environment.R;  // the residual vectors

        if (R.cols() > 0) {  // if there are residual vectors available
            const auto are_any_values_larger = (R.colwise().norm().array() > this->threshold).any();
            return !are_any_values_larger;
        }

        return false;  // there aren't any residual vectors available
    }
};


}  // namespace GQCP
