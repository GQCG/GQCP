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
#include "QCMethod/Geminals/GeminalCoefficientsInterface.hpp"


namespace GQCP {


/*
 *  PUBLIC METHODS
 */

/**
 *  @param fock_space       the seniority-zero ONV basis the wave function should live in
 *
 *  @return the wave function expansion corresponding to the geminal coefficients
 */
LinearExpansion GeminalCoefficientsInterface::toLinearExpansion(const SpinUnresolvedONVBasis& fock_space) const {

    // The ONV can't be marked const, as makeONV() and setNext() are non-const methods

    VectorX<double> coefficients = VectorX<double>::Zero(fock_space.get_dimension());  // coefficient vector
    SpinUnresolvedONV onv = fock_space.makeONV(0);  // start with address 0
    for (size_t I = 0; I < fock_space.get_dimension(); I++) {

        coefficients(I) = this->overlap(onv);

        if (I < fock_space.get_dimension() - 1) {  // skip the last permutation
            fock_space.setNextONV(onv);
        }
    }

    return LinearExpansion(fock_space, coefficients);
}


}  // namespace GQCP
