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
#include "Processing/Properties/BaseElectricalResponseSolver.hpp"

#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationSolver.hpp"


namespace GQCP {



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the linear response equations for the wave function response
 * 
 *  @param sq_hamiltonian           the Hamiltonian parameters expressed in an orthonormal orbital basis
 *  @param dipole_op                the dipole integrals expressed in an orthonormal orbital basis
 * 
 *  @return the wave function response as an (Nx3)-matrix
 */
Matrix<double, Dynamic, 3> BaseElectricalResponseSolver::calculateWaveFunctionResponse(const SQHamiltonian<double>& sq_hamiltonian, const VectorSQOneElectronOperator<double> dipole_op) const {

    const auto k_p = this->calculateParameterResponseConstant(sq_hamiltonian);  // p for parameter
    const auto F_p = this->calculateParameterResponseForce(dipole_op);  // has 3 columns


    // This function is basically a wrapper around solving k_p x = -F_p
    auto environment = LinearEquationEnvironment<double>(k_p, -F_p);
    auto solver = LinearEquationSolver<double>::HouseholderQR();
    solver.perform(environment);

    const Matrix<double, Dynamic, 3> x = environment.x;
    return x;
}


}  // namespace GQCP
