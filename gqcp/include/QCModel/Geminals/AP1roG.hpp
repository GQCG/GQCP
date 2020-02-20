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


/**
 *  The AP1roG (=pCCD) geminal (coupled-cluster) wave function model.
 */
class AP1roG {

public:

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
     *  @param G                the AP1roG geminal coefficients
     *  @param multipliers      the AP1roG Lagrangian multipliers
     *
     *  @return the AP1roG response 1-DM
     */
    static OneRDM<double> calculate1RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers);

    /**
     *  @param G                the AP1roG geminal coefficients
     *  @param multipliers      the AP1roG Lagrangian multipliers
     *
     *  @return the AP1roG response number 2-RDM (the Delta-matrix in the notes)
     */
    static SquareMatrix<double> calculateNumber2RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers);

    /**
     *  @param G                the AP1roG geminal coefficients
     *  @param multipliers      the AP1roG Lagrangian multipliers
     *
     *  @return the AP1roG response pair 2-RDM (the Pi-matrix in the notes)
     */
    static SquareMatrix<double> calculatePair2RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers);

    /**
     *  @param G                the AP1roG geminal coefficients
     *  @param multipliers      the AP1roG Lagrangian multipliers
     *
     *  @return the AP1roG response 2-DM
     */
    static TwoRDM<double> calculate2RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers);
};



}  // namespace GQCP
