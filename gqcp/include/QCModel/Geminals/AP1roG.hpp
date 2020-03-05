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


#include "Mathematical/Representation/BlockMatrix.hpp"
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
        G (G)
    {}


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
};



}  // namespace QCModel
}  // namespace GQCP
