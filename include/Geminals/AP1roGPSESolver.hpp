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


#include "Geminals/BaseAP1roGSolver.hpp"

#include "Geminals/AP1roGPSEs.hpp"


namespace GQCP {


/**
 *  A class that is able to solve the AP1roG projected Schrödinger equations (PSEs)
 */
class AP1roGPSESolver {
    double convergence_threshold;
    size_t maximum_number_of_iterations;

    AP1roGPSEs pses;  // the AP1roG PSEs


public:
    // CONSTRUCTORS

    /**
     *  @param pses         the AP1roG PSEs
     */
    AP1roGPSESolver(const AP1roGPSEs& pses, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128);


    // PUBLIC METHODS

    /**
     *  Solve the projected Schrödinger equations for AP1roG
     * 
     *  @param G            the initial geminal coefficients, that are updated in every iteration to the converged geminal coefficients
     */
    void solve(AP1roGGeminalCoefficients& G) const;

    /**
     *  Solve the projected Schrödinger equations for AP1roG, using a zero initial guess
     */
    AP1roGGeminalCoefficients solve() const;
};


}  // namespace GQCP
