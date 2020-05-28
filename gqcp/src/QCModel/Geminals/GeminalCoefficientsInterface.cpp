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

#include "QCModel/Geminals/GeminalCoefficientsInterface.hpp"


namespace GQCP {


/*
 *  PUBLIC METHODS
 */

/**
 *  @param onv_basis       the seniority-zero spin-resolved ONV basis the wave function should live in
 *
 *  @return the wave function expansion corresponding to the geminal coefficients
 */
LinearExpansion<SeniorityZeroONVBasis> GeminalCoefficientsInterface::toLinearExpansion(const SeniorityZeroONVBasis& onv_basis) const {

    const auto dim = onv_basis.dimension();
    VectorX<double> coefficients = VectorX<double>::Zero(dim);  // the expansion coefficient vector

    const auto onv_basis_proxy = onv_basis.proxy();
    SpinUnresolvedONV onv = onv_basis_proxy.constructONVFromAddress(0);  // start with address 0
    for (size_t I = 0; I < dim; I++) {

        coefficients(I) = this->overlap(onv);

        if (I < dim - 1) {  // prevent the last permutation from occurring
            onv_basis_proxy.transformONVToNextPermutation(onv);
        }
    }

    return LinearExpansion<SeniorityZeroONVBasis>(onv_basis, coefficients);
}


}  // namespace GQCP
