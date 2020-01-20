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


#include "Mathematical/Optimization/Algorithm/IterationCycle.hpp"


#include <cstddef>


namespace GQCP {



/**
 *  An algorithm that performs iterations. In every iteration, convergence is checked and a set of iteration steps is performed.
 * 
 *  @tparam _Iterate            the type of the iterate
 */
template <typename Environment>
class IterativeAlgorithm {
public:

private:
    size_t maximum_number_of_iterations;
    size_t iteration = 0;  // the current iteration counter

    IterationCycle<Environment> iteration_cycle;
    


public:


    /*
     *  CONSTRUCTORS
     */
    
    /**
     *  Initialize the members of the iterative algorithm
     * 
     *  @param maximum_number_of_iterations         the maximum number of iterations the algorithm may perform
     */
    IterativeSolver(const size_t maximum_number_of_iterations = 128) :
        maximum_number_of_iterations (maximum_number_of_iterations),
    {}


    /**
     *  Perform the iteration steps 
     */
    void iterate(Environment& environment) {

        for (this->iteration = 0; this->iteration < this->maximum_number_of_iterations; this->iteration++) {



        }

    }




};



}  // namespace GQCP
