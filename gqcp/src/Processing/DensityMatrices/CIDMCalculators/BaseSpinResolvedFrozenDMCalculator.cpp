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

#include "Processing/DensityMatrices/CIDMCalculators/BaseSpinResolvedFrozenDMCalculator.hpp"

#include "Utilities/linalg.hpp"


namespace GQCP {

/*
 *  CONSTRUCTORS
 */

/**
 *  @param dm_calculator                shared pointer to active (non-frozen core) DM builder
 *  @param X                            the number of frozen orbitals
 */
BaseSpinResolvedFrozenDMCalculator::BaseSpinResolvedFrozenDMCalculator(const std::shared_ptr<BaseSpinResolvedDMCalculator> dm_calculator, const size_t X) :
    BaseSpinResolvedDMCalculator(),
    active_dm_calculator {std::move(dm_calculator)},
    X {X} {}


/*
 * PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 1-DMs given a coefficient vector
 */
SpinResolvedOneDM<double> BaseSpinResolvedFrozenDMCalculator::calculateSpinResolved1DM(const VectorX<double>& x) const {

    auto K = this->onvBasis()->numberOfOrbitals();

    OneDM<double> D_aa = OneDM<double>::Zero(K, K);
    OneDM<double> D_bb = OneDM<double>::Zero(K, K);

    auto K_active = K - this->X;

    // Set the values for the frozen orbital
    for (size_t i = 0; i < this->X; i++) {
        D_aa(i, i) = 1;
        D_bb(i, i) = 1;
    }

    SpinResolvedOneDM<double> sub_1DMs = this->active_dm_calculator->calculateSpinResolved1DM(x);

    // Incorporate the submatrices from the active space
    D_aa.block(this->X, this->X, K_active, K_active) += sub_1DMs.alpha();
    D_bb.block(this->X, this->X, K_active, K_active) += sub_1DMs.beta();

    return SpinResolvedOneDM<double>(D_aa, D_bb);
};


/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 2-DMs given a coefficient vector
 */
SpinResolvedTwoDM<double> BaseSpinResolvedFrozenDMCalculator::calculateSpinResolved2DM(const VectorX<double>& x) const {

    auto K = this->onvBasis()->numberOfOrbitals();

    TwoDM<double> d_aaaa {K};
    d_aaaa.setZero();

    TwoDM<double> d_aabb {K};
    d_aabb.setZero();

    TwoDM<double> d_bbaa {K};
    d_bbaa.setZero();

    TwoDM<double> d_bbbb {K};
    d_bbbb.setZero();


    SpinResolvedOneDM<double> one_DMs = this->active_dm_calculator->calculateSpinResolved1DM(x);
    auto D_aa = one_DMs.alpha();
    auto D_bb = one_DMs.beta();

    // Implement frozen DM formulas
    for (size_t p = 0; p < this->X; p++) {  // iterate over frozen orbitals

        // DM Overlap between frozen and active space:
        //      frozen orbital indices (p) must always have one annihilation and one creation index (always occupied)
        //      values are dictated by the 'active' orbital indices and correspond to that of the active 1DMs
        //      Hence we start adding the 1DMs starting from index 'X' the number frozen orbitals
        d_aaaa.addBlock<0, 1>(D_aa, this->X, this->X, p, p);
        d_aaaa.addBlock<2, 3>(D_aa, p, p, this->X, this->X);
        d_aaaa.addBlock<2, 1>(-D_aa, p, this->X, this->X, p);
        d_aaaa.addBlock<3, 0>(-D_aa, this->X, p, p, this->X);

        d_bbbb.addBlock<0, 1>(D_bb, this->X, this->X, p, p);
        d_bbbb.addBlock<2, 3>(D_bb, p, p, this->X, this->X);
        d_bbbb.addBlock<2, 1>(-D_bb, p, this->X, this->X, p);
        d_bbbb.addBlock<3, 0>(-D_bb, this->X, p, p, this->X);

        d_aabb.addBlock<2, 3>(D_bb, p, p, this->X, this->X);
        d_aabb.addBlock<0, 1>(D_aa, this->X, this->X, p, p);

        d_bbaa.addBlock<2, 3>(D_aa, p, p, this->X, this->X);
        d_bbaa.addBlock<0, 1>(D_bb, this->X, this->X, p, p);


        // Set the values for the frozen orbital
        d_bbaa(p, p, p, p) = 1;
        d_aabb(p, p, p, p) = 1;

        for (size_t q = p + 1; q < this->X; q++) {  // iterate over frozen orbitals

            d_aaaa(p, p, q, q) = 1;
            d_aaaa(q, q, p, p) = 1;
            d_aaaa(p, q, q, p) = -1;
            d_aaaa(q, p, p, q) = -1;

            d_bbbb(p, p, q, q) = 1;
            d_bbbb(q, q, p, p) = 1;
            d_bbbb(p, q, q, p) = -1;
            d_bbbb(q, p, p, q) = -1;

            d_aabb(p, p, q, q) = 1;
            d_bbaa(p, p, q, q) = 1;

            d_aabb(q, q, p, p) = 1;
            d_bbaa(q, q, p, p) = 1;
        }
    }


    // Incorporate the 2-DM subblocks into the total 2DMs
    SpinResolvedTwoDM<double> sub_2DMs = this->active_dm_calculator->calculateSpinResolved2DM(x);

    d_aaaa.addBlock(sub_2DMs.alphaAlpha(), this->X, this->X, this->X, this->X);
    d_bbbb.addBlock(sub_2DMs.betaBeta(), this->X, this->X, this->X, this->X);
    d_aabb.addBlock(sub_2DMs.alphaBeta(), this->X, this->X, this->X, this->X);
    d_bbaa.addBlock(sub_2DMs.betaAlpha(), this->X, this->X, this->X, this->X);

    return SpinResolvedTwoDM<double>(d_aaaa, d_aabb, d_bbaa, d_bbbb);
};


}  // namespace GQCP
