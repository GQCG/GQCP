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


#include "Mathematical/Optimization/IterativeSolver.hpp"


#include <type_traits>


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
class AcceleratedIterativeSolver :
    public IterativeSolver<_Iterate> {

public:
    using Iterate = _Iterate;
    using Accelerator = _Accelerator;
    using Proxy = typename Accelerator::Subject;
    using BaseSolver = IterativeSolver<_Iterate>;


protected:
    Accelerator accelerator;  // the accelerator used for the acceleration of the iterates


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize the solver with an initial guess and an accelerator
     * 
     *  @param initial_guess                        the initial guess to the solver
     *  @param accelerator                          the accelerator used for the acceleration of the iterates
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform
     */
    AcceleratedIterativeSolver(const Iterate& initial_guess, const Accelerator& accelerator, const size_t maximum_number_of_iterations = 128) :
        BaseSolver(initial_guess, maximum_number_of_iterations),
        accelerator (accelerator)
    {}


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
     *  Calculate the next iterate as if there was no acceleration on it
     * 
     *  @return the 'regular' iterate
     */
    virtual Iterate calculateRegularIterate() = 0;

    /**
     *  Calculate a next iterate from its proxy
     * 
     * @return the next iterated, calculated from its proxy
     */
    virtual Iterate calculateIterateFromProxy(const Proxy& proxy) = 0;


    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @return a new iterate to be used in the next iteration
     */
    Iterate calculateNextIterate() {

        // Calculate a next iterate through accelerating its proxy
        const auto iterate = this->calculateRegularIterate();
        const auto accelerator_subject = this->calculateAcceleratorSubject(iterate);
        const auto accelerated_subject = this->accelerator.accelerate(accelerator_subject);

        return this->calculateIterateFromProxy(accelerated_subject);
    }
};


}  // namespace GQCP
