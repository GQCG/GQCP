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


/**
 *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
 *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
 *
 *  @return an element of the N-DM, as specified by the given bra and ket indices
 *
 *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
 */
double GeneralSpinUnresolvedDMCalculator::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) const {

    if (this->coefficients.rows() == 0) {
        throw std::logic_error("GeneralSpinUnresolvedDMCalculator::calculateElement(std::vector<size_t>, std::vector<size_t>): No vector has been set.");
    }

    return this->dm_calculator.calculateElement(bra_indices, ket_indices, this->coefficients);
}


}  // namespace GQCP
