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


#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"


namespace GQCP {


/**
 *  A class that represents a wave function: expansion coefficients in a spin-unresolved ONV basis
 */
class SpinUnresolvedLinearExpansion : public LinearExpansion {

    /**
     *  @param onv_basis            the spin-unresolved ONV basis in which the wave function 'lives'
     *  @param coefficients         the expansion coefficients
     */
    SpinUnresolvedLinearExpansion(const SpinUnresolvedONVBasis& onv_basis, const VectorX<double>& coefficients);
};


}  // namespace GQCP
