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
 *  A protocol (virtual base class) for quantum chemical methods
 * 
 *  @tparam InitialGuess            the type of the initial guess that is supplied to the solution algorithm
 *  @tparam Solution                the type of the solution that is returned by the solution algorithm
 */
template <typename InitialGuess, typename Solution>
class QCMethodProtocol {

    /**
     *  Determine the optimal wave function parameters for the quantum chemical method
     * 
     *  @return the optimal wave function parameters. This is often encapsulated with other information, like the energy
     */
    virtual Solution solve(const InitialGuess& initial_guess) = 0;
};



}  // namespace GQCP
