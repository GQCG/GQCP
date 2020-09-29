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


#include "DensityMatrix/CIDMCalculators/BaseSpinResolvedFrozenDMCalculator.hpp"
#include "DensityMatrix/SpinResolvedOneDM.hpp"
#include "DensityMatrix/SpinResolvedTwoDM.hpp"
#include "ONVBasis/SpinResolvedFrozenONVBasis.hpp"


namespace GQCP {


/**
 *  A class capable of calculating 1- and 2-DMs from wave functions expanded in the frozen DOCI ONV basis
 */
class FrozenCoreDOCIRDMBuilder: public BaseSpinResolvedFrozenDMCalculator {
private:
    SpinUnresolvedFrozenONVBasis onv_basis;  // both the frozen alpha and beta spin-unresolved ONV basis


public:
    // CONSTRUCTORS

    /**
     *  @param onv_basis        both the frozen alpha and beta spin-unresolved ONV basis
     */
    explicit FrozenCoreDOCIRDMBuilder(const SpinUnresolvedFrozenONVBasis& onv_basis);


    // OVERRIDDEN PUBLIC METHODS

    /**
     *  @return the ONV basis that is associated to this RDMBuilder
     */
    const BaseONVBasis* onvBasis() const override { return &onv_basis; }
};


}  // namespace GQCP
