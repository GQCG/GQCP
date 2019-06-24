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
#ifndef GQCP_AP1ROGPSES_HPP
#define GQCP_AP1ROGPSES_HPP


#include "Geminals/AP1roGGeminalCoefficients.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Mathematical/BlockMatrix.hpp"
#include "Mathematical/BlockRankFourTensor.hpp"
#include "typedefs.hpp"


namespace GQCP {


/**
 *  A class that represents the AP1roG projected Schrödinger equations (PSEs)
 */
class AP1roGPSEs {
private:
    size_t N_P;  // the number of electron pairs
    size_t K;  // the number of spatial orbitals

    HamiltonianParameters<double> ham_par;  // the one- and two-electron integrals in an orthonormal orbital basis


public:
    // CONSTRUCTORS

    /**
     *  @param ham_par          the one- and two-electron integrals in an orthonormal orbital basis
     *  @param N_P              the number of electron pairs
     */
    AP1roGPSEs(const HamiltonianParameters<double>& ham_par, const size_t N_P);


    // PUBLIC METHODS

    /**
     *  @return the number of electron pairs
     */
    size_t numberOfElectronPairs() const { return this->N_P; }

    /**
     *  @return the number of electron pairs
     */
    size_t numberOfSpatialOrbitals() const { return this->K; }

    /**
     *  @param G        the AP1roG geminal coefficients
     *  @param i        the subscript for the coordinate function
     *  @param a        the superscript for the coordinate function
     *
     *  @return the value of the coordinate function with given indices (i,a) at the given geminal coefficients
     */
    double evaluateCoordinateFunction(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a) const;

    /**
     *  @param G            the geminal coefficients
     * 
     *  @return the PSEs, evaluated at the given geminal coefficients
     */
    BlockMatrix<double> evaluateCoordinateFunctions(const AP1roGGeminalCoefficients& G) const;

    /**
     *  @return a callable (i.e. with operator()) expression for the coordinate functions
     */
    VectorFunction callableCoordinateFunctions() const;

    /**
     *  @param G        the AP1roG geminal coefficients
     *  @param i        the subscript for the coordinate function
     *  @param a        the superscript for the coordinate function
     *  @param j        the subscript for the geminal coefficient
     *  @param b        the superscript for the geminal coefficient
     *
     *  @return the value of the Jacobian element with compound indices (i,a) and (j,b) at the given geminal coefficients
     */
    double evaluateJacobianElement(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a, const size_t j, const size_t b) const;

    /**
     *  @param G            the geminal coefficients
     * 
     *  @return the Jacobian, J_{ia,jb} of the PSEs (df_i^a/dG_j^b), evaluated at the given geminal coefficients
     */
    BlockRankFourTensor<double> evaluateJacobian(const AP1roGGeminalCoefficients& G) const;

    /**
     *  @return a callable expression for the Jacobian
     */
    MatrixFunction callableJacobian() const;
};


}  // namespace GQCP



#endif  // GQCP_AP1ROGPSES_HPP
