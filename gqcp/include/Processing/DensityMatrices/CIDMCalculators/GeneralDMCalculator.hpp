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


#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/BaseSpinResolvedDMCalculator.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/adaptor/strided.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>

#include <memory>


namespace GQCP {


/**
 *  A wrapper around the derived DMCalculators that provides the functionality of the appropriate derived DMCalculator for a given ONV basis at compile- or runtime.
 */
class GeneralDMCalculator {
private:
    std::shared_ptr<BaseSpinResolvedDMCalculator> dm_calculator;
    VectorX<double> coefficients;

public:
    // CONSTRUCTOR

    /**
     *  The default constructor.
     */
    GeneralDMCalculator() = default;

    /**
     *  Allocate an SpinResolvedDMCalculator
     *
     *  @param onv_basis       the FCI ONV basis
     */
    explicit GeneralDMCalculator(const SpinResolvedONVBasis& onv_basis);

    /**
     *  Allocate a SpinResolvedSelectedDMCalculator
     *
     *  @param onv_basis       the 'selected' ONV basis
     */
    explicit GeneralDMCalculator(const SpinResolvedSelectedONVBasis& onv_basis);

    /**
     *  A run-time constructor allocating the appropriate derived DMCalculator
     *
     *  @param onv_basis       the ONV basis on which the DMCalculator should be based
     */
    explicit GeneralDMCalculator(const BaseONVBasis& onv_basis);

    /**
     *  A run-time constructor allocating the appropriate derived DMCalculator and coefficient vector
     *
     *  @param linear_expansion       the linear expansion in a certain ONV basis
     */
    template <typename ONVBasis>
    explicit GeneralDMCalculator(const LinearExpansion<ONVBasis>& linear_expansion) :
        GeneralDMCalculator(linear_expansion.onvBasis()) {

        this->setCoefficients(linear_expansion.coefficients());
    }


    // PUBLIC METHODS

    /**
     *  @return all 1-DMs if a given coefficient vector is set
     */
    SpinResolvedOneDM<double> calculateSpinResolved1DM() const;

    /**
     *  @return all 2-DMs if a given coefficient vector is set
     */
    SpinResolvedTwoDM<double> calculateSpinResolved2DM() const;

    /**
     *  Replace this instance's coefficients with the given coefficients.
     * 
     *  @param coefficients                 the new coefficients for this instance
     */
    void setCoefficients(const VectorX<double>& coefficients) { this->coefficients = coefficients; };
};


}  // namespace GQCP
