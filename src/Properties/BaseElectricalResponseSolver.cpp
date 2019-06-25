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
#include "Properties/BaseElectricalResponseSolver.hpp"

#include "Eigen/Householder"


namespace GQCP {



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the linear response equations for the wave function response
 * 
 *  @param ham_par                  the Hamiltonian parameters
 *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
 * 
 *  @return the wave function response
 */
Matrix<double, Dynamic, 3> BaseElectricalResponseSolver::calculateWaveFunctionResponse(const HamiltonianParameters<double>& ham_par, const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const {

    const auto k_p = this->calculateParameterResponseConstant(ham_par);  // p for parameter
    const auto F_p = this->calculateParameterResponseForce(dipole_integrals);  // has 3 columns


    // This function is basically a wrapper around solving k_p x = -F_p
    Eigen::HouseholderQR<Eigen::MatrixXd> linear_solver (k_p);
    const Matrix<double, Dynamic, 3> x = linear_solver.solve(-F_p);

    return x;
}


}  // namespace GQCP
