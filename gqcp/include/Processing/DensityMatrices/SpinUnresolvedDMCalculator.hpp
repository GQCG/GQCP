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
#include "Processing/DensityMatrices/BaseSpinUnresolvedDMCalculator.hpp"
#include "Processing/DensityMatrices/SpinResolvedOneDM.hpp"
#include "Processing/DensityMatrices/SpinResolvedTwoDM.hpp"


namespace GQCP {


/**
 *  A class capable of calculating RDMs from wave functions expanded in the full CI spin-unresolved ONV basis
 */
class SpinUnresolvedDMCalculator: public BaseSpinUnresolvedDMCalculator {
private:
    SpinUnresolvedONVBasis onv_basis;  // spin-unresolved ONV basis


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param onv_basis                spin-unresolved ONV basis
     */
    explicit SpinUnresolvedDMCalculator(const SpinUnresolvedONVBasis& onv_basis);

    /**
     *  The default constructor.
     */
    SpinUnresolvedDMCalculator() = default;


    // DESTRUCTOR

    /**
     *  The default destructor.
     */
    ~SpinUnresolvedDMCalculator() = default;


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param x        the coefficient vector representing the UnresolvedCI wave function
     *
     *  @return the 1-RDM given a coefficient vector
     */
    OneDM<double> calculate1RDM(const VectorX<double>& x) const override { throw std::runtime_error("SpinUnresolvedDMCalculator::calculate1RDMs(VectorX<double>): not implemented yet"); }

    /**
     *  @param x        the coefficient vector representing the UnresolvedCI wave function
     *
     *  @return the 2-RDM given a coefficient vector
     */
    TwoDM<double> calculate2RDM(const VectorX<double>& x) const override { throw std::runtime_error("SpinUnresolvedDMCalculator::calculate2RDMs(VectorX<double>): not implemented yet"); }

    /**
     *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
     *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
     *  @param x                the coefficient vector representing the UnresolvedCI wave function
     *
     *  @return an element of the N-RDM, as specified by the given bra and ket indices
     *
     *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
     */
    double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const VectorX<double>& x) const override;

    /**
     *  @return the ONV basis that is associated to this RDMBuilder
     */
    const BaseONVBasis* onvBasis() const override { return &onv_basis; }
};


}  // namespace GQCP
