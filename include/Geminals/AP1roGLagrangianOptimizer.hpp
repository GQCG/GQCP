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


#include "Geminals/AP1roGGeminalCoefficients.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

namespace GQCP {


/**
 *  A class that is able to determine the Lagrange multipliers for the AP1roG PSE Lagrangian, given a solution for the geminal coefficients
 */
class AP1roGLagrangianOptimizer {
private:
    AP1roGGeminalCoefficients G;  // the converged geminal coefficients that are a solution to the AP1roG PSEs
    SQHamiltonian<double> sq_hamiltonian;  // the Hamiltonian expressed in an orthonormal orbital basis


public:

    // CONSTRUCTORS

    /**
     *  @param G                    the converged geminal coefficients that are a solution to the AP1roG PSEs
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal orbital basis
     */
    AP1roGLagrangianOptimizer(const AP1roGGeminalCoefficients& G, const SQHamiltonian<double>& sq_hamiltonian);


    // PUBLIC METHODS

    /**
     *  @return the Lagrange multipliers for the AP1roG PSE Lagrangian
     */
    AP1roGVariables solve();
};


}  // namespace GQCP
