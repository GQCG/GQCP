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
#include "QCModel/Geminals/AP1roG.hpp"


namespace GQCP {


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param G                    the converged AP1roG geminal coefficients
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the AP1roG electronic energy
 */
double QCModel::AP1roG::calculateEnergy(const AP1roGGeminalCoefficients& G, const SQHamiltonian<double>& sq_hamiltonian) {

    // Prepare some variables
    const auto& h = sq_hamiltonian.core().parameters();
    const auto& g = sq_hamiltonian.twoElectron().parameters();


    // KISS implementation of the AP1roG energy
    double E = 0.0;
    for (size_t j = 0; j < G.get_N_P(); j++) {
        E += 2 * h(j,j);

        for (size_t k = 0; k < G.get_N_P(); k++) {
            E += 2 * g(k,k,j,j) - g(k,j,j,k);
        }

        for (size_t b = G.get_N_P(); b < G.get_K(); b++) {
            E += g(j,b,j,b) * G(j,b);
        }
    }

    return E;
}


}  // namespace GQCP
