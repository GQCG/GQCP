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
#ifndef GQCP_AP1ROGPSESOLVER_HPP
#define GQCP_AP1ROGPSESOLVER_HPP


#include "Geminals/AP1roGPSEs.hpp"


namespace GQCP {


/**
 *  A class that is able to solve the AP1roG projected Schrödinger equations (PSEs)
 */
class AP1roGPSESolver {
    AP1roGPSEs pses;  // the AP1roG PSEs


public:
    // CONSTRUCTORS

    /**
     *  @param pses         the AP1roG PSEs
     */
    AP1roGPSESolver(const AP1roGPSEs& pses);


    // PUBLIC METHODS

    /**
     *  Solve the projected Schrödinger equations for AP1roG
     * 
     *  @param G            the initial geminal coefficients
     */
    AP1roGGeminalCoefficients solve(const AP1roGGeminalCoefficients& G) const;

    /**
     *  Solve the projected Schrödinger equations for AP1roG, using a zero initial guess
     */
    AP1roGGeminalCoefficients solve() const;
};


}  // namespace GQCP



#endif  // GQCP_AP1ROGPSESOLVER_HPP
