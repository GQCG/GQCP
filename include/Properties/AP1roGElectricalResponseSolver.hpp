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
#ifndef GQCP_AP1ROGELECTRICALRESPONSESOLVER_HPP
#define GQCP_AP1ROGELECTRICALRESPONSESOLVER_HPP


#include "Properties/BaseElectricalResponseSolver.hpp"


namespace GQCP {


class AP1roGElectricalResponseSolver : public BaseElectricalResponseSolver {
public:
    // CONSTRUCTORS


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param ham_par                  the Hamiltonian parameters
     * 
     *  @return the parameter response constant (k_p), i.e. the first-order parameter partial derivative of the PSEs function (i.e. the Jacobian of the PSEs)
     */
    SquareMatrix<double> calculateParameterResponseConstant(const HamiltonianParameters<double>& ham_par) const override;

    /**
     *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
     * 
     *  @return the parameter response force (F_p), i.e. the first-order perturbation derivative of the PSEs
     */
    Matrix<double, Dynamic, 3> calculateParameterResponseForce(const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const override;


    // PUBLIC METHODS

    /**
     *  @param ham_par                  the Hamiltonian parameters
     * 
     *  @return the Lagrangian multiplier response constant (k_lambda), i.e. the first-order perturbation derivative of the PSEs
     */
    SquareMatrix<double> calculateMultiplierResponseConstant(const HamiltonianParameters<double>& ham_par) const;

    /**
     *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
     * 
     *  @return the Lagrangian multiplier response force (F_p)
     */
    Matrix<double, Dynamic, 3> calculateMultiplierResponseForce(const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const;

    /**
     *  Solve the linear response equations for the Lagrangian multiplier response
     * 
     *  @param ham_par                  the Hamiltonian parameters
     *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
     * 
     *  @return the Lagrangian multiplier
     */
    Matrix<double, Dynamic, 3> calculateMultiplierResponse(const HamiltonianParameters<double>& ham_par, const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const;
};


}  // namespace GQCP



#endif  // GQCP_AP1ROGELECTRICALRESPONSESOLVER_HPP
