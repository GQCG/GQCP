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


#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/BaseSpinUnresolvedDMCalculator.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SpinUnresolvedDMCalculator.hpp"

#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/adaptor/strided.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>

#include <memory>


namespace GQCP {


/**
 *  A wrapper around the SpinUnresolvedDMCalculator
 */
class GeneralSpinUnresolvedDMCalculator {
private:
    SpinUnresolvedDMCalculator dm_calculator;
    VectorX<double> coefficients;


public:
    // CONSTRUCTORS

    /**
     *  The default constructor.
     */
    GeneralSpinUnresolvedDMCalculator() = default;

    /**
     *  Allocate a SpinUnresolvedDMCalculator
     *
     *  @param onv_basis       the spin-unresolved ONV basis
     */
    explicit GeneralSpinUnresolvedDMCalculator(const SpinUnresolvedONVBasis& onv_basis);


    // PUBLIC METHODS

    /**
     *  @return the 1-DM if a given coefficient vector is set
     */
    OneDM<double> calculate1DM() const;

    /**
     *  @return the 2-DM if a given coefficient vector is set
     */
    TwoDM<double> calculate2DM() const;

    /**
     *  Replace this instance's coefficients with new ones.
     * 
     *  @param coefficients                 the new coefficients
     */
    void setCoefficients(const VectorX<double>& coefficients) { this->coefficients = coefficients; };
};


}  // namespace GQCP
