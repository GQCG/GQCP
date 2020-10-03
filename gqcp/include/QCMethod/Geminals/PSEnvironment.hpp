// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationEnvironment.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/Geminals/AP1roG.hpp"


namespace GQCP {
namespace PSEnvironment {


/**
 *  Create a non-linear equation environment with the AP1roG geminal coefficients as an initial guess.
 * 
 *  @param sq_hamiltonian                   the second-quantized Hamiltonian, expressed in an orthonormal spinor basis
 *  @param G_initial                        the initial guess for the AP1roG geminal coefficients
 * 
 *  @return an environment suitable for AP1roG calculations
 */
template <typename Scalar>
NonLinearEquationEnvironment<Scalar> AP1roG(const RSQHamiltonian<Scalar>& sq_hamiltonian, const AP1roGGeminalCoefficients& G_initial) {

    const auto N_P = G_initial.numberOfElectronPairs();

    const auto initial_guess = G_initial.asVector();  // column major
    const auto f_callable = QCModel::AP1roG::callablePSECoordinateFunctions(sq_hamiltonian, N_P);
    const auto J_callable = QCModel::AP1roG::callablePSEJacobian(sq_hamiltonian, N_P);

    return GQCP::NonLinearEquationEnvironment<Scalar>(initial_guess, f_callable, J_callable);
}


/**
 *  Create a non-linear equation environment with a zero initial guess for the AP1roG geminal coefficients.
 * 
 *  @param sq_hamiltonian                   the second-quantized Hamiltonian, expressed in an orthonormal spinor basis
 *  @param N_P                              the number of electron pairs
 * 
 *  @return an environment suitable for AP1roG calculations
 */
template <typename Scalar>
NonLinearEquationEnvironment<Scalar> AP1roG(const RSQHamiltonian<Scalar>& sq_hamiltonian, const size_t N_P) {

    const auto K = sq_hamiltonian.numberOfOrbitals();    // number of spatial orbitals
    const AP1roGGeminalCoefficients G_initial {N_P, K};  // geminal coefficients set to zero

    return GQCP::PSEnvironment::AP1roG<Scalar>(sq_hamiltonian, G_initial);
}


}  // namespace PSEnvironment
}  // namespace GQCP
