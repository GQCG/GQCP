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
#include "DensityMatrix/SpinResolved1DM.hpp"
#include "DensityMatrix/SpinResolvedTwoDM.hpp"
#include "ONVBasis/SpinResolvedFrozenONVBasis.hpp"


namespace GQCP {


/**
 *  A class capable of calculating 1- and 2-DMs from wave functions expanded in the full frozen full spin resolved ONV basis
 */
class SpinResolvedFrozenDMCalculator: public BaseSpinResolvedFrozenDMCalculator {
private:
    SpinResolvedFrozenONVBasis onv_basis;


public:
    // CONSTRUCTORS

    /**
     *  @param onv_basis       the frozen spin-resolved ONV basis
     */
    explicit SpinResolvedFrozenDMCalculator(const SpinResolvedFrozenONVBasis& onv_basis);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @return the ONV basis that is associated to this DMCalculator
     */
    const BaseONVBasis* onvBasis() const override { return &onv_basis; }
};


}  // namespace GQCP
