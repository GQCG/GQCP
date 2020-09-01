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


#include "ONVBasis/BaseONVBasis.hpp"
#include "Processing/DensityMatrices/OneDM.hpp"
#include "Processing/DensityMatrices/SpinResolvedOneDM.hpp"
#include "Processing/DensityMatrices/SpinResolvedTwoDM.hpp"
#include "Processing/DensityMatrices/TwoDM.hpp"


namespace GQCP {


/**
 *  BaseSpinResolvedDMCalculator is an abstract base class whose derived classes are capable of calculating 1- and 2-DMs from wave functions in a spin resolved basis
 */
class BaseSpinResolvedDMCalculator {
public:
    // CONSTRUCTORS
    BaseSpinResolvedDMCalculator() = default;


    // DESTRUCTOR
    virtual ~BaseSpinResolvedDMCalculator() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 1-DMs given a coefficient vector
     */
    virtual SpinResolvedOneDM<double> calculateSpinResolved1DM(const VectorX<double>& x) const = 0;

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 2-DMs given a coefficient vector
     */
    virtual SpinResolvedTwoDM<double> calculateSpinResolved2DM(const VectorX<double>& x) const = 0;

    /**
     *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
     *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
     *  @param x                the coefficient vector representing the wave function
     *
     *  @return an element of the spin-summed (total) N-DM, as specified by the given bra and ket indices.
     *
     *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
     */
    virtual double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const VectorX<double>& x) const = 0;

    /**
     *  @return the ONV basis associated to this DMCalculator
     */
    virtual const BaseONVBasis* onvBasis() const = 0;
};


}  // namespace GQCP
