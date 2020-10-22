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


#include "Processing/Properties/BaseElectricalResponseSolver.hpp"


namespace GQCP {


/**
 *  A class that is able to solve the electrical CPHF-equations for RHF
 */
class RHFElectricalResponseSolver: public BaseElectricalResponseSolver {
private:
    size_t N_P;  // the number of electron pairs


public:
    // CONSTRUCTORS

    /**
     *  @param N_P          the number of electron pairs
     */
    RHFElectricalResponseSolver(const size_t N_P);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param sq_hamiltonian           the Hamiltonian parameters
     * 
     *  @return the parameter response constant (k_p), i.e. the second-order parameter partial derivative of the RHF energy function
     */
    SquareMatrix<double> calculateParameterResponseConstant(const RSQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  @param dipole_op                the dipole integrals expressed in an orthonormal orbital basis
     * 
     *  @return the parameter response force (F_p) as an (Nx3)-matrix, i.e. the first-order parameter partial derivative of the perturbation derivative of the RHF energy function
     */
    Matrix<double, Dynamic, 3> calculateParameterResponseForce(const VectorRSQOneElectronOperator<double>& dipole_op) const override;
};


}  // namespace GQCP
