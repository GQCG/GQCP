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

#include "Processing/DensityMatrices/CIDMCalculators/GeneralDMCalculator.hpp"

#include "Processing/DensityMatrices/CIDMCalculators/SeniorityZeroDMCalculator.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SpinResolvedDMCalculator.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SpinResolvedSelectedDMCalculator.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */


/**
 *  Allocate a SpinResolvedDMCalculator
 *
 *  @param onv_basis       the FCI ONV basis
 */
GeneralDMCalculator::GeneralDMCalculator(const SpinResolvedONVBasis& onv_basis) :
    dm_calculator {std::make_shared<SpinResolvedDMCalculator>(onv_basis)} {}


/**
 *  Allocate a SpinResolvedSelectedDMCalculator
 *
 *  @param onv_basis       the 'selected' ONV basis
 */
GeneralDMCalculator::GeneralDMCalculator(const SpinResolvedSelectedONVBasis& onv_basis) :
    dm_calculator {std::make_shared<SpinResolvedSelectedDMCalculator>(onv_basis)} {}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return all 1-DMs if a given coefficient vector is set
 */
SpinResolvedOneDM<double> GeneralDMCalculator::calculateSpinResolved1DM() const {

    if (this->coefficients.rows() == 0) {
        throw std::logic_error("GeneralDMCalculator::calculateSpinResolved1DM(): No vector has been set.");
    }

    return dm_calculator->calculateSpinResolved1DM(this->coefficients);
}


/**
 *  @return all 2-DMs if a given coefficient vector is set
 */
SpinResolvedTwoDM<double> GeneralDMCalculator::calculateSpinResolved2DM() const {

    if (this->coefficients.rows() == 0) {
        throw std::logic_error("GeneralDMCalculator::calculateSpinResolved2DM(): No vector has been set.");
    }

    return dm_calculator->calculateSpinResolved2DM(this->coefficients);
}


}  // namespace GQCP
