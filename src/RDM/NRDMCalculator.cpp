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
#include "RDM/NRDMCalculator.hpp"



#include <iostream>



namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the FCI Fock space with spin orbitals
 *  @param coeff            the expansion coefficient vector
 */
NRDMCalculator::NRDMCalculator(const FockSpace& fock_space, const Eigen::VectorXd& coeff) :
    fock_space (fock_space),
    coeff (coeff)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
 *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
 *
 *  @return an element of the N-RDM, as specified by the given bra and ket indices
 *
 *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2)
 */
double NRDMCalculator::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) const {


    // The ket indices should be reversed because the annihilators on the ket should be applied from right to left
    std::vector<size_t> ket_indices_reversed = ket_indices;
    std::reverse(ket_indices_reversed.begin(), ket_indices_reversed.end());


    double value = 0.0;
    int bra_sign = 1;
    int ket_sign = 1;

    auto it_bra = fock_space.begin();
    auto end = fock_space.end();  // help the compiler by putting the end out of the for-loop
    while (it_bra != end) {

        bra_sign = 1;
        ONV bra = it_bra.currentONV();

        // Annihilate the bra on the bra indices
        if (!bra.annihilateAll(bra_indices, bra_sign)) {
            ++it_bra;
            continue;
        }

        auto it_ket = fock_space.begin();
        while (it_ket != end) {

            ket_sign = 1;
            ONV ket = it_ket.currentONV();

            // Annihilate the ket on the ket indices
            if (!ket.annihilateAll(ket_indices_reversed, ket_sign)) {
                ++it_ket;
                continue;
            }

            size_t I = it_bra.currentAddress();
            size_t J = it_ket.currentAddress();
            if (bra == ket) {
                value += bra_sign * ket_sign * this->coeff(I) * this->coeff(J);
            }

            ++it_ket;
        }  // ket Fock space iteration

        ++it_bra;
    }  // bra Fock space iteration

    return value;
}



}  // namespace GQCP
