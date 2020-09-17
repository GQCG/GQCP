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

#include "Processing/DensityMatrices/CIDMCalculators/GeneralSpinUnresolvedDMCalculator.hpp"

#include "Processing/DensityMatrices/CIDMCalculators/SeniorityZeroDMCalculator.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SpinResolvedDMCalculator.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SpinResolvedSelectedDMCalculator.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SpinUnresolvedDMCalculator.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Allocate a SpinUnresolvedDMCalculator
 *
 *  @param onv_basis       the ONV basis
 */
GeneralSpinUnresolvedDMCalculator::GeneralSpinUnresolvedDMCalculator(const SpinUnresolvedONVBasis& onv_basis) :
    dm_calculator {SpinUnresolvedDMCalculator(onv_basis)} {}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the 1-DM if a given coefficient vector is set
 */
OneDM<double> GeneralSpinUnresolvedDMCalculator::calculate1DM() const {

    if (this->coefficients.rows() == 0) {
        throw std::logic_error("GeneralSpinUnresolvedDMCalculator::calculate1DM(): No vector has been set.");
    }

    return dm_calculator.calculate1DM(this->coefficients);
}

/**
 *  @return the 2-DM if a given coefficient vector is set
 */
TwoDM<double> GeneralSpinUnresolvedDMCalculator::calculate2DM() const {

    if (this->coefficients.rows() == 0) {
        throw std::logic_error("GeneralSpinUnresolvedDMCalculator::calculate2DM(): No vector has been set.");
    }

    return dm_calculator.calculate2DM(this->coefficients);
}


}  // namespace GQCP
