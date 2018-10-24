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
#include "AP1roG/AP1roG.hpp"


namespace GQCP {

/*
 *  CONSTRUCTORS
 */
/**
 *  Default constructor setting everything to zero
 */
AP1roG::AP1roG() :
    geminal_coefficients (AP1roGGeminalCoefficients()),
    electronic_energy (0.0)
{}


/**
 *  Constructor based on given @param geminal_coefficients and @param electronic_energy
 */
AP1roG::AP1roG(const GQCP::AP1roGGeminalCoefficients& geminal_coefficients, double electronic_energy) :
    geminal_coefficients (geminal_coefficients),
    electronic_energy (electronic_energy)
{}



/*
 *  HELPER FUNCTIONS
 */
/**
 *  Calculate the AP1roG energy given AP1roG geminal coefficients @param G and Hamiltonian parameters @param ham_par
 */
double calculateAP1roGEnergy(const GQCP::AP1roGGeminalCoefficients& G, const GQCP::HamiltonianParameters& ham_par) {

    Eigen::MatrixXd h_SO = ham_par.get_h().get_matrix_representation();
    Eigen::Tensor<double, 4> g_SO = ham_par.get_g().get_matrix_representation();


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


}  // namespace GQCP
