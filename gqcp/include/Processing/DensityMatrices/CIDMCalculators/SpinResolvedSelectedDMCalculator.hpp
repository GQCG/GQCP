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


#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "Processing/DensityMatrices/SpinResolvedOneDM.hpp"
#include "Processing/DensityMatrices/SpinResolvedTwoDM.hpp"


namespace GQCP {


/**
 *  A class capable of calculating 1- and 2-DMs from wave functions expanded in a selected spin-resolved basis
 */
class SpinResolvedSelectedDMCalculator {
private:
    SpinResolvedSelectedONVBasis onv_basis;  // spin-resolved ONV basis containing the selected configurations


public:
    // CONSTRUCTORS

    /**
     *  @param onv_basis                spin-resolved ONV basis containing the selected configurations
     */
    explicit SpinResolvedSelectedDMCalculator(const SpinResolvedSelectedONVBasis& onv_basis);


    // DESTRUCTOR

    /**
     *  The default destructor.
     */
    ~SpinResolvedSelectedDMCalculator() = default;


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param x        the coefficient vector representing the 'selected' wave function
     *
     *  @return all 1-DMs given a coefficient vector
     */
    SpinResolvedOneDM<double> calculateSpinResolved1DM(const VectorX<double>& x) const;

    /**
     *  @param x        the coefficient vector representing the 'selected' wave function
     *
     *  @return all 2-DMs given a coefficient vector
     */
    SpinResolvedTwoDM<double> calculateSpinResolved2DM(const VectorX<double>& x) const;
};


}  // namespace GQCP
