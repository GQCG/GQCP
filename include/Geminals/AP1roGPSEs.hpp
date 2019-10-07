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
 *  A class that represents the actual AP1roG PSEs: it contains formulas for the evaluation of the PSEs at the given AP1roG geminal coefficients, as well as the Jacobian at the given AP1roG geminal coefficients
 */
class AP1roGPSEs {
private:
    size_t N_P;  // the number of electron pairs
    size_t K;  // the number of spatial orbitals

    SQHamiltonian<double> sq_hamiltonian;  // the one- and two-electron integrals in an orthonormal orbital basis


public:
    // CONSTRUCTORS

    /**
     *  @param sq_hamiltonian       the one- and two-electron integrals in an orthonormal orbital basis
     *  @param N_P                  the number of electron pairs
     */
    AP1roGPSEs(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P);


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
     *  @return the coordinate function with given indices (i,a) at the given geminal coefficients
     */
    double calculateCoordinateFunction(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a) const;

    /**
     *  @param G        the AP1roG geminal coefficients
     *
     *  @return a row-major vector of coordinate functions at the given geminal coefficients
     */
    VectorX<double> calculateCoordinateFunctions(const AP1roGGeminalCoefficients& G) const;

    /**
     *  @return a callable (i.e. with operator()) expression for the coordinate functions
     */
    VectorFunction callableCoordinateFunctions() const;

    /**
     *  @param G        the AP1roG geminal coefficients
     *  @param i        the subscript for the coordinate function
     *  @param a        the superscript for the coordinate function
     *  @param k        the subscript for the geminal coefficient
     *  @param c        the superscript for the geminal coefficient
     *
     *  @return the Jacobian element with compound indices (i,a) and (k,c) at the given geminal coefficients
     */
    double calculateJacobianElement(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a, const size_t k, const size_t c) const;

    /**
     *  @param G        the AP1roG geminal coefficients
     *
     *  @return the Jacobian (in a row-major representation) at the given geminal coefficients
     */
    SquareMatrix<double> calculateJacobian(const AP1roGGeminalCoefficients& G) const;

    /**
     *  @return a callable expression for the Jacobian
     */
    MatrixFunction callableJacobian() const;
};


}  // namespace GQCP
