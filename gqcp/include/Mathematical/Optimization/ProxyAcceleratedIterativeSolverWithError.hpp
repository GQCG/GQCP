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


#include "Mathematical/Optimization/ProxyAcceleratedIterativeSolver.hpp"


namespace GQCP {


/**
 *  An iterative solver that uses an accelerator on proxies of the iterates to help convergence. It partially implements calculateNextIterate(), but requires derived classes to implement:
 *      - calculateRegularIterate();
 *      - calculateAcceleratorSubject();
 *      - calculateIterateFromProxy();
 * 
 *  @tparam _Iterate            the type of the iterate
 *  @tparam _Accelerator        the type of the accelerator that is used
 */
template <typename _Iterate, typename _Accelerator>
class ProxyAcceleratedIterativeSolver :
    public BaseAcceleratedIterativeSolver<_Iterate> {


public:
    using Iterate = _Iterate;
    using Accelerator = _Accelerator;
    using Error = typename Accelerator::Error;
    using Proxy = typename Accelerator::Subject;
    using Base = BaseAcceleratedIterativeSolver<Iterate>;


public:

    /*
     *  CONSTRUCTORS
     */

    using BaseAcceleratedIterativeSolver<Iterate>::BaseAcceleratedIterativeSolver;  // inherit base constructors


    /*
     *  PURE VIRTUAL PUBLIC METHODS
     */

    /**
     *  Calculate the proxy that corresponds to the iterate, which can be used by the accelerator
     * 
     *  @return the proxy that corresponds to the iterate
     */
    virtual Proxy calculateAcceleratorSubject(const Iterate& iterate) = 0;

    /**
     *  Calculate a next iterate from its proxy
     * 
     * @return the next iterated, calculated from its proxy
     */
    virtual Iterate calculateIterateFromProxy(const Proxy& proxy) = 0;


    /**
     *  @return the error measure of the current iteration
     */
    virtual Error calculateError() = 0;


    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @return a new iterate to be used in the next iteration
     */
    Iterate calculateNextIterate() {

        // Calculate a next iterate through accelerating its proxy and corresponding error
        const auto iterate = this->calculateRegularIterate();
        const auto accelerator_subject = this->calculateAcceleratorSubject(iterate);
        const auto error = this->calculateError();
        this->accelerator.feed(accelerator_subject, error);
        const auto accelerated_subject = this->accelerator.accelerate();

        return this->calculateIterateFromProxy(accelerated_subject);
    }
};


}  // namespace GQCP
