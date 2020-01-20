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


#include "Mathematical/Optimization/BaseAcceleratedIterativeSolver.hpp"

#include <type_traits>


namespace GQCP {


/**
 *  An iterative solver that uses an accelerator on the iterates to help convergence. It partially implements calculateNextIterate(), but requires derived classes to implement calculateRegularIterate().
 * 
 *  @tparam _Iterate            the type of the iterate
 *  @tparam _Accelerator        the type of the accelerator that is used
 */
template <typename _Iterate, typename _Accelerator>
class AcceleratedIterativeSolver :
    public BaseAcceleratedIterativeSolver<_Iterate> {
static_assert(std::is_same<_Iterate, typename _Accelerator::Subject>);

public:
    using Iterate = _Iterate;
    using Accelerator = _Accelerator;
    using Base = BaseAcceleratedIterativeSolver<Iterate>;


public:

    /*
     *  CONSTRUCTORS
     */

    using BaseAcceleratedIterativeSolver<Iterate>::BaseAcceleratedIterativeSolver;  // inherit base constructors


    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @return a new iterate to be used in the next iteration
     */
    virtual Iterate calculateNextIterate() {

        const auto iterate = this->calculateRegularIterate();
        this->accelerator.feed(iterate);
        return this->accelerator.accelerate(iterate);
    }
};


}  // namespace GQCP
