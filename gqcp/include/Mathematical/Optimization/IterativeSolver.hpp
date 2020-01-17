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
 *  An iterative solver: at every step, the iterate is updated
 */
template <typename _Iterate>
class IterativeSolver {
private:
    Iterate iterate;  // the iterate, which hopefully converges to the solution


public:

    /*
     *  CONSTRUCTORS
     */
    
    /**
     *  Initialize the solver with an initial guess
     * 
     *  @param initial_guess            the initial guess to the solver
     */
    IterativeSolver(const Iterate& initial_guess) : 
        iterate (initial_guess)
    {}



    /*
     *  PUBLIC METHODS
     */
    
    /**
     *  Iterate until the algorithm has converged
     * 
     *  @return the solution of the algorithm
     */
    Iterate solve() {
        while (!this->derived().isConverged()) {
            this->iterate = this->derived().updateIterate();
        }
    }
};


}  // namespace GQCP
