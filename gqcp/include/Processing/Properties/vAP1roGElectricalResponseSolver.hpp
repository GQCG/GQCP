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

#include "Processing/Properties/BaseElectricalResponseSolver.hpp"

#include "QCModel/Geminals/vAP1roG.hpp"


namespace GQCP {


/**
 *  A class whose instances can solve the response equations for vAP1roG.
 */
class vAP1roGElectricalResponseSolver : public BaseElectricalResponseSolver {
private:
    QCModel::vAP1roG vap1rog;  // the optimal vAP1roG parameters


public:
    // CONSTRUCTORS

    /**
     *  @param vap1rog          the optimal vAP1roG parameters
     */
    vAP1roGElectricalResponseSolver(const QCModel::vAP1roG& vap1rog);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param sq_hamiltonian                  the Hamiltonian expressed in an orthonormal orbital basis
     * 
     *  @return the parameter response constant (k_p), i.e. the first-order parameter partial derivative of the PSEs, which is the Jacobian of the PSEs
     */
    SquareMatrix<double> calculateParameterResponseConstant(const SQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  @param dipole_op                the dipole integrals expressed in an orthonormal orbital basis
     * 
     *  @return the parameter response force (F_p) as an (Nx3)-matrix, i.e. the first-order perturbation derivative of the PSEs
     */
    Matrix<double, Dynamic, 3> calculateParameterResponseForce(const VectorSQOneElectronOperator<double>& dipole_op) const override;


    // PUBLIC METHODS

    /**
     *  @param sq_hamiltonian                   the Hamiltonian expressed in an orthonormal orbital basis
     * 
     *  @return the Lagrangian multiplier response constant (k_lambda), which is the transpose of the parameter multiplier response constant
     */
    SquareMatrix<double> calculateMultiplierResponseConstant(const SQHamiltonian<double>& sq_hamiltonian) const;

    /**
     *  @param dipole_op                        the dipole integrals expressed in an orthonormal orbital basis
     * 
     *  @return the explicit (i.e. the first part of the) Lagrangian multiplier response force, A_lambda
     */
    Matrix<double, Dynamic, 3> calculateExplicitMultiplierResponseForce(const VectorSQOneElectronOperator<double> dipole_op) const;

    /**
     *  @param sq_hamiltonian                   the Hamiltonian expressed in an orthonormal orbital basis
     *
     *  @return the multiplier force constant of the implicit part (i.e. the second part of the) Lagrangian multiplier response, B_lambda
     */
    BlockRankFourTensor<double> calculateImplicitMultiplierResponseForceConstant(const SQHamiltonian<double>& sq_hamiltonian) const;

    /**
     *  @param sq_hamiltonian                   the Hamiltonian expressed in an orthonormal orbital basis
     *  @param dipole_op                        the dipole integrals expressed in an orthonormal orbital basis
     *  @param x                                the first-order parameter response
     * 
     *  @return the Lagrangian multiplier response force (F_lambda)
     */
    Matrix<double, Dynamic, 3> calculateMultiplierResponseForce(const SQHamiltonian<double>& sq_hamiltonian, const VectorSQOneElectronOperator<double> dipole_op, const Matrix<double, Dynamic, 3>& x) const;

    /**
     *  Solve the linear response equations for the Lagrangian multiplier response
     * 
     *  @param sq_hamiltonian               the Hamiltonian parameters
     *  @param dipole_op                    the dipole integrals expressed in an orthonormal orbital basis
     *  @param x                            the first-order parameter response
     * 
     *  @return the Lagrangian multiplier reponse y
     */
    Matrix<double, Dynamic, 3> calculateMultiplierResponse(const SQHamiltonian<double>& sq_hamiltonian, const VectorSQOneElectronOperator<double> dipole_op, const Matrix<double, Dynamic, 3>& x) const;
};


}  // namespace GQCP
