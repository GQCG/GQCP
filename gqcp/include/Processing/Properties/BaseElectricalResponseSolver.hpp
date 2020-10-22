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


#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  A class that serves as a base class to solve electrical linear response equations
 */
class BaseElectricalResponseSolver {
public:
    // DESTRUCTOR
    virtual ~BaseElectricalResponseSolver() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  @param sq_hamiltonian           the Hamiltonian in an orthonormal orbital basis
     * 
     *  @return the parameter response constant (k_p)
     */
    virtual SquareMatrix<double> calculateParameterResponseConstant(const RSQHamiltonian<double>& sq_hamiltonian) const = 0;

    /**
     *  @param dipole_op                the dipole integrals expressed in an orthonormal orbital basis
     * 
     *  @return the parameter response force (F_p) as an (Nx3)-matrix
     */
    virtual Matrix<double, Dynamic, 3> calculateParameterResponseForce(const VectorRSQOneElectronOperator<double>& dipole_op) const = 0;


    // PUBLIC METHODS

    /**
     *  Solve the linear response equations for the wave function response
     * 
     *  @param sq_hamiltonian           the Hamiltonian parameters expressed in an orthonormal orbital basis
     *  @param dipole_op                the dipole integrals expressed in an orthonormal orbital basis
     * 
     *  @return the wave function response as an (Nx3)-matrix
     */
    Matrix<double, Dynamic, 3> calculateWaveFunctionResponse(const RSQHamiltonian<double>& sq_hamiltonian, const VectorRSQOneElectronOperator<double> dipole_op) const;
};


}  // namespace GQCP
