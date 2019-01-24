// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#include "geminals/AP1roG.hpp"


namespace GQCP {


/**
 *  @param G            the converged AP1roG geminal coefficients
 *  @param ham_par      Hamiltonian parameters in an orthonormal spatial orbital basis
 *
 *  @return the AP1roG electronic energy
 */
double calculateAP1roGEnergy(const AP1roGGeminalCoefficients& G, const HamiltonianParameters& ham_par) {

    OneElectronOperator h_SO = ham_par.get_h();
    TwoElectronOperator g_SO = ham_par.get_g();


    // KISS implementation of the AP1roG energy
    double E = 0.0;
    for (size_t j = 0; j < G.get_N_P(); j++) {
        E += 2 * h_SO(j,j);

        for (size_t k = 0; k < G.get_N_P(); k++) {
            E += 2 * g_SO(k,k,j,j) - g_SO(k,j,j,k);
        }

        for (size_t b = G.get_N_P(); b < G.get_K(); b++) {
            E += g_SO(j,b,j,b) * G(j,b);
        }
    }

    return E;
}


/**
 *  @param G            the AP1roG geminal coefficients
 *  @param Q            the AP1roG bivariational coefficients
 *
 *  @return the overlap between the bivariational coefficients and the geminal coefficients, i.e. <Phi(q)|Psi(p)>
 */
double calculateOverlap(const AP1roGGeminalCoefficients& G, const BivariationalCoefficients& Q) {

    double overlap = Q.q0;
    AP1roGVariables q = Q.q;

    for (size_t i = 0; i < G.get_N_P(); i++) {
        for (size_t a = G.get_N_P(); a < G.get_K(); a++) {
            overlap += q(i,a) * G(i,a);
        }
    }

    return overlap;
}



}  // namespace GQCP
