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
#ifndef GQCP_BASEELECTRICALRESPONSESOLVER_HPP
#define GQCP_BASEELECTRICALRESPONSESOLVER_HPP


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Mathematical/Matrix.hpp"
#include "Operator/OneElectronOperator.hpp"

#include <array>



namespace GQCP {


class BaseElectricalResponseSolver {
public:
    // DESTRUCTOR
    virtual ~BaseElectricalResponseSolver() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  @param ham_par                  the Hamiltonian parameters
     * 
     *  @return the parameter response constant (k_p)
     */
    virtual SquareMatrix<double> calculateParameterResponseConstant(const HamiltonianParameters<double>& ham_par) const = 0;

    /**
     *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
     * 
     *  @return the parameter response force (F_p)
     */
    virtual Matrix<double, Dynamic, 3> calculateParameterResponseForce(const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const = 0;


    // PUBLIC METHODS

    /**
     *  Solve the linear response equations for the wave function response
     * 
     *  @param ham_par                  the Hamiltonian parameters
     *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
     * 
     *  @return the wave function response
     */
    Matrix<double, Dynamic, 3> calculateWaveFunctionResponse(const HamiltonianParameters<double>& ham_par, const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const;
};


}  // namespace GQCP



#endif  // GQCP_BASEELECTRICALRESPONSESOLVER_HPP
