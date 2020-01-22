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
 *  A standard environment that stores the iterative variables of an iterative algorithm
 * 
 *  @tparam _Iterate            the type of the iterative variables
 */
template <typename _Iterate>
class OptimizationEnvironment {
public:
    using Iterate = _Iterate;
    using Scalar = typename Iterate::Scalar;


public:
    std::deque<_Iterate> variables;  // a collection of variables that iteratively grows through an optimization algorithm; for example all iterates x when solving f(x) = 0 iteratively


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize the optimization environment with an initial guess
     * 
     *  @param initial_guess                the initial guess for the variables
     */
    OptimizationEnvironment(const Iterate& initial_guess) : 
        variables(1, initial_guess)
    {}
};


}  // namespace GQCP
