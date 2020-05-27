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


#include "Mathematical/Representation/ImplicitMatrixSlice.hpp"
#include "Mathematical/Representation/ImplicitRankFourTensorSlice.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/Geminals/AP1roGGeminalCoefficients.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The AP1roG (=pCCD) geminal (coupled-cluster) wave function model.
 */
class AP1roG {
private:
    AP1roGGeminalCoefficients G;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param G                the optimal geminal coefficients
     */
    AP1roG(const AP1roGGeminalCoefficients& G) :
        G {G} {}


    /*
     *  STATIC PUBLIC METHODS
     */

    /**
     *  @param G                    the AP1roG geminal coefficients
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the AP1roG electronic energy
     */
    static double calculateEnergy(const AP1roGGeminalCoefficients& G, const SQHamiltonian<double>& sq_hamiltonian);

    /**
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param G                        the AP1roG geminal coefficients
     *  @param i                        the subscript for the coordinate function
     *  @param a                        the superscript for the coordinate function
     *
     *  @return the PSE coordinate function with given indices (i,a) at the given geminal coefficients
     */
    static double calculatePSECoordinateFunction(const SQHamiltonian<double>& sq_hamiltonian, const AP1roGGeminalCoefficients& G, const size_t i, const size_t a);

    /**
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param G                        the AP1roG geminal coefficients
     *
     *  @return the PSEs, evaluated at the given geminal coefficients
     */
    static ImplicitMatrixSlice<double> calculatePSECoordinateFunctions(const SQHamiltonian<double>& sq_hamiltonian, const AP1roGGeminalCoefficients& G);

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param G                    the AP1roG geminal coefficients
     *
     *  @return the Jacobian J_{ia,jb} of the PSEs, i.e. df_i^a/dG_j^b, evaluated at the given geminal coefficients
     */
    static ImplicitRankFourTensorSlice<double> calculatePSEJacobian(const SQHamiltonian<double>& sq_hamiltonian, const AP1roGGeminalCoefficients& G);

    /**
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param N_P                      the number of electron pairs
     * 
     *  @return a callable (i.e. with operator()) expression for the Jacobian: the accepted VectorX<double> argument should contain the geminal coefficients in a column-major representation
     */
    static MatrixFunction<double> callablePSEJacobian(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P);

    /**
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param N_P                      the number of electron pairs
     * 
     *  @return a callable (i.e. with operator()) expression for the coordinate functions: the accepted VectorX<double> argument should contain the geminal coefficients in a column-major representation
     */
    static VectorFunction<double> callablePSECoordinateFunctions(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P);

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param G                    the AP1roG geminal coefficients
     *  @param i                    the subscript for the coordinate function
     *  @param a                    the superscript for the coordinate function
     *  @param j                    the subscript for the geminal coefficient
     *  @param b                    the superscript for the geminal coefficient
     *
     *  @return the Jacobian element with compound indices (i,a) and (j,b) at the given geminal coefficients
     */
    static double calculatePSEJacobianElement(const SQHamiltonian<double>& sq_hamiltonian, const AP1roGGeminalCoefficients& G, const size_t i, const size_t a, const size_t j, const size_t b);


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     * 
     *  @return the electronic energy for these AP1roG model parameters
     */
    double calculateEnergy(const SQHamiltonian<double>& sq_hamiltonian) const { return AP1roG::calculateEnergy(this->G, sq_hamiltonian); }

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the Jacobian J_{ia,jb} of the PSEs, i.e. df_i^a/dG_j^b, evaluated at these AP1roG model parameters
     */
    ImplicitRankFourTensorSlice<double> calculatePSEJacobian(const SQHamiltonian<double>& sq_hamiltonian) const { return AP1roG::calculatePSEJacobian(sq_hamiltonian, this->G); }

    /**
     *  @return the corresponding geminal coefficients of these AP1roG model parameters
     */
    const AP1roGGeminalCoefficients& geminalCoefficients() const { return this->G; }

    /**
     *  @return the number of electron pairs that are described by these AP1roG model parameters
     */
    size_t numberOfElectronPairs() const { return this->G.numberOfElectronPairs(); }

    /**
     *  @return the number of spatial orbitals that are described by these AP1roG model parameters
     */
    size_t numberOfSpatialOrbitals() const { return this->G.numberOfSpatialOrbitals(); }

    /**
     *  @return the implicit occupied-virtual orbital space that is associated with these AP1roG model parameters
     */
    OrbitalSpace orbitalSpace() const { return this->geminalCoefficients().orbitalSpace(); }
};


}  // namespace QCModel
}  // namespace GQCP
