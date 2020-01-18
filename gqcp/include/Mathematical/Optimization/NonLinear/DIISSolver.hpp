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

#include <deque>


namespace GQCP {


/**
 *  A base solver that accelerates its iterates using the direct inversion of the iterative subspace.
 * 
 *  @tparam _Iterate                the type of the iterate
 *  @tparam 
 *  @tparam _Error                  the type of the error estimate
 * 
 *  This class achieves compile-time polymorphism through CRTP. Derived classes should implement:
 *      - isConverged(), to check if the algorithm is considered to be converged
 *      - regularNewIterate(), to produce a new iterate to be used in the next iteration if the DIIS acceleration has not been switched on yet
 */
template <typename _Iterate, typename _, typename _Error>
class DIISSolver : public IterativeSolver<_Iterate> {
public:
    using Iterate = _Iterate;
    using Base = IterativeSolver<_Iterate>;


private:
    size_t minimum_subspace_dimension;  // the minimum number of iterates that have to be in the subspace before enabling the DIIS acceleration
    size_t maximum_subspace_dimension;  // the maximum DIIS subspace dimension before the oldest Fock matrices get discarded (one at a time)



public:

    /**
     *  @param minimum_subspace_dimension       the minimum number of iterates that have to be in the subspace before enabling the DIIS acceleration
     *  @param maximum_subspace_dimension       the maximum DIIS subspace dimension before the oldest Fock matrices get discarded (one at a time)
     */
    DIISSolver(const size_t minimum_subspace_dimension, const size_t maximum_subspace_dimension) :
        minimum_subspace_dimension (minimum_subspace_dimension),
        maximum_subspace_dimension (maximum_subspace_dimension)
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     * 
     */
    Iterate calculateDIISIterate() {

        // Start by creating a regular new iterate
        Iterate iterate = this->regularNewIterate();


        // Calculate a better iterate (using DIIS acceleration)

    }


    /**
     *  @return a new iterate according to the DIIS algorithm
     */
    Iterate updateIterate() {

        if (this->subspaceDimension() < this->minimum_subspace_dimension) {
            // Don't do DIIS
            return this->regularNewIterate();
        }

        else {
            // Enable the DIIS acceleration
            return this->calculateDIISIterate();
        }
    }
};


}  // namespace GQCP
