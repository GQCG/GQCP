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
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/Geminals/AP1roG.hpp"
#include "QCModel/Geminals/AP1roGGeminalCoefficients.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The AP1roG (=pCCD) geminal (coupled-cluster) wave function model that is variationally determined
 */
class vAP1roG {
private:
    AP1roGGeminalCoefficients G;              // the optimized geminal coefficients
    ImplicitMatrixSlice<double> multipliers;  // the optimized Lagrange multipliers


public:
    // CONSTRUCTORS

    /**
     *  @param G                    the optimal geminal coefficients
     *  @param multipliers          the optimized Lagrange multipliers
     */
    vAP1roG(const AP1roGGeminalCoefficients& G, const ImplicitMatrixSlice<double>& multipliers);


    //  STATIC PUBLIC METHODS

    /**
     *  @param G                the AP1roG geminal coefficients
     *  @param multipliers      the AP1roG Lagrangian multipliers
     *
     *  @return the AP1roG response 1-DM
     */
    static OneDM<double> calculate1RDM(const AP1roGGeminalCoefficients& G, const ImplicitMatrixSlice<double>& multipliers);

    /**
     *  @param G                the AP1roG geminal coefficients
     *  @param multipliers      the AP1roG Lagrangian multipliers
     *
     *  @return the AP1roG response 2-DM
     */
    static TwoDM<double> calculate2RDM(const AP1roGGeminalCoefficients& G, const ImplicitMatrixSlice<double>& multipliers);

    /**
     *  @param G                    the AP1roG geminal coefficients
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the AP1roG electronic energy
     */
    static double calculateEnergy(const AP1roGGeminalCoefficients& G, const SQHamiltonian<double>& sq_hamiltonian) { return AP1roG::calculateEnergy(G, sq_hamiltonian); }

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param N_P                  the number of electron pairs
     * 
     *  @return the response force (-F_lambda) that is used to solve the linear equations for the Lagrange multipliers lambda in [k_lambda lambda = -F_lambda]
     */
    static ImplicitMatrixSlice<double> calculateMultiplierResponseForce(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P);

    /**
     *  @param G                    the AP1roG geminal coefficients
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     * 
     *  @return the response force constant (k_lambda) that is used to solve the linear equations for the Lagrange multipliers lambda in [k_lambda lambda = -F_lambda]
     */
    static MatrixX<double> calculateMultiplierResponseForceConstant(const SQHamiltonian<double>& sq_hamiltonian, const AP1roGGeminalCoefficients& G);

    /**
     *  @param G                the AP1roG geminal coefficients
     *  @param multipliers      the AP1roG Lagrangian multipliers
     *
     *  @return the AP1roG response number 2-RDM (the Delta-matrix in the notes)
     */
    static SquareMatrix<double> calculateNumber2RDM(const AP1roGGeminalCoefficients& G, const ImplicitMatrixSlice<double>& multipliers);

    /**
     *  @param G                the AP1roG geminal coefficients
     *  @param multipliers      the AP1roG Lagrangian multipliers
     *
     *  @return the AP1roG response pair 2-RDM (the Pi-matrix in the notes)
     */
    static SquareMatrix<double> calculatePair2RDM(const AP1roGGeminalCoefficients& G, const ImplicitMatrixSlice<double>& multipliers);


    // PUBLIC METHODS

    /**
     *  @return the reponse one-electron density matrix for these vAP1roG parameters
     */
    OneDM<double> calculate1RDM() const { return vAP1roG::calculate1RDM(this->G, this->multipliers); };

    /**
     *  @return the reponse two-electron density matrix for these vAP1roG parameters
     */
    TwoDM<double> calculate2RDM() const { return vAP1roG::calculate2RDM(this->G, this->multipliers); };

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     * 
     *  @return the electronic energy for these AP1roG model parameters
     */
    double calculateEnergy(const SQHamiltonian<double>& sq_hamiltonian) const { return vAP1roG::calculateEnergy(this->G, sq_hamiltonian); }

    /**
     *  @return the corresponding geminal coefficients of these AP1roG model parameters
     */
    const AP1roGGeminalCoefficients& geminalCoefficients() const { return this->G; }

    /**
     *  @return the Lagrange multipliers
     */
    const ImplicitMatrixSlice<double>& lagrangeMultipliers() const { return this->multipliers; }

    /**
     *  @return the number of electron pairs that are described by these AP1roG model parameters
     */
    size_t numberOfElectronPairs() const { return this->G.numberOfElectronPairs(); }

    /**
     *  @return the number of spatial orbitals that are described by these AP1roG model parameters
     */
    size_t numberOfSpatialOrbitals() const { return this->G.numberOfSpatialOrbitals(); }
};


}  // namespace QCModel
}  // namespace GQCP
