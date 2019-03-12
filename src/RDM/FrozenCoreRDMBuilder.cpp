// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "RDM/FrozenCoreRDMBuilder.hpp"
#include "utilities/linalg.hpp"


namespace GQCP {

/*
 *  CONSTRUCTORS
 */

/**
 *  @param rdm_builder                  shared pointer to active (non-frozen core) RDM builder
 *  @param X                            the number of frozen orbitals
 */
FrozenCoreRDMBuilder::FrozenCoreRDMBuilder(std::shared_ptr<BaseRDMBuilder> rdm_builder, size_t X) :
    BaseRDMBuilder(),
    active_rdm_builder (std::move(rdm_builder)),
    X (X)
{}



/*
 * OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 1-RDMs given a coefficient vector
 */
OneRDMs<double> FrozenCoreRDMBuilder::calculate1RDMs(const VectorX<double>& x) const {

    auto K = this->get_fock_space()->get_K();

    OneRDM<double> D_aa = OneRDM<double>::Zero(K, K);
    OneRDM<double> D_bb = OneRDM<double>::Zero(K, K);

    auto K_active = K - this->X;

    // Set the values for the frozen orbital
    for (size_t i = 0; i < this->X; i++) {
        D_aa(i,i) = 1;
        D_bb(i,i) = 1;
    }

    OneRDMs<double> sub_1rdms = this->active_rdm_builder->calculate1RDMs(x);

    // Incorporate the submatrices from the active space
    D_aa.block(this->X, this->X, K_active, K_active) += sub_1rdms.one_rdm_aa;
    D_bb.block(this->X, this->X, K_active, K_active) += sub_1rdms.one_rdm_bb;

    return OneRDMs<double>(D_aa, D_bb);
};


/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 2-RDMs given a coefficient vector
 */
TwoRDMs<double> FrozenCoreRDMBuilder::calculate2RDMs(const VectorX<double>& x) const {

    auto K = this->get_fock_space()->get_K();

    TwoRDM<double> d_aaaa (K);
    d_aaaa.setZero();

    TwoRDM<double> d_aabb (K);
    d_aabb.setZero();

    TwoRDM<double> d_bbaa (K);
    d_bbaa.setZero();

    TwoRDM<double> d_bbbb (K);
    d_bbbb.setZero();


    OneRDMs<double> one_rdms = this->active_rdm_builder->calculate1RDMs(x);
    auto D_aa = one_rdms.one_rdm_aa;
    auto D_bb = one_rdms.one_rdm_bb;

    // Implement frozen RDM formulas
    for (size_t p = 0; p < this->X; p++) {  // iterate over frozen orbitals

        // RDM Overlap between frozen and active space:
        //      frozen orbital indices (p) must always have one annihilation and one creation index (always occupied)
        //      values are dictated by the 'active' orbital indices and correspond to that of the active 1RDMs
        //      Hence we start adding the 1RDMs starting from index 'X' the number frozen orbitals
        d_aaaa.addBlock<0,1>(D_aa, this->X, this->X, p, p);
        d_aaaa.addBlock<2,3>(D_aa, p, p, this->X, this->X);
        d_aaaa.addBlock<2,1>(-D_aa, p, this->X, this->X, p);
        d_aaaa.addBlock<3,0>(-D_aa, this->X, p, p, this->X);

        d_bbbb.addBlock<0,1>(D_bb, this->X, this->X, p, p);
        d_bbbb.addBlock<2,3>(D_bb, p, p, this->X, this->X);
        d_bbbb.addBlock<2,1>(-D_bb, p, this->X, this->X, p);
        d_bbbb.addBlock<3,0>(-D_bb, this->X, p, p, this->X);

        d_aabb.addBlock<2,3>(D_bb, p, p, this->X, this->X);
        d_aabb.addBlock<0,1>(D_aa, this->X, this->X, p, p);

        d_bbaa.addBlock<2,3>(D_aa, p, p, this->X, this->X);
        d_bbaa.addBlock<0,1>(D_bb, this->X, this->X, p, p);


        // Set the values for the frozen orbital
        d_bbaa(p,p,p,p) = 1;
        d_aabb(p,p,p,p) = 1;

        for (size_t q = p +1; q < this->X; q++) {  // iterate over frozen orbitals

            d_aaaa(p,p,q,q) = 1;
            d_aaaa(q,q,p,p) = 1;
            d_aaaa(p,q,q,p) = -1;
            d_aaaa(q,p,p,q) = -1;

            d_bbbb(p,p,q,q) = 1;
            d_bbbb(q,q,p,p) = 1;
            d_bbbb(p,q,q,p) = -1;
            d_bbbb(q,p,p,q) = -1;

            d_aabb(p,p,q,q) = 1;
            d_bbaa(p,p,q,q) = 1;

            d_aabb(q,q,p,p) = 1;
            d_bbaa(q,q,p,p) = 1;
        }
    }


    // Incorporate the 2-RDM subblocks into the total 2RDMs
    TwoRDMs<double> sub_2rdms = this->active_rdm_builder->calculate2RDMs(x);

    d_aaaa.addBlock(sub_2rdms.two_rdm_aaaa, this->X, this->X, this->X, this->X);
    d_bbbb.addBlock(sub_2rdms.two_rdm_bbbb, this->X, this->X, this->X, this->X);
    d_aabb.addBlock(sub_2rdms.two_rdm_aabb, this->X, this->X, this->X, this->X);
    d_bbaa.addBlock(sub_2rdms.two_rdm_bbaa, this->X, this->X, this->X, this->X);

    return TwoRDMs<double>(d_aaaa, d_aabb, d_bbaa, d_bbbb);
};


/**
 *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
 *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
 *  @param x                the coefficient vector representing the wave function
 *
 *  @return an element of the spin-summed (total) N-RDM, as specified by the given bra and ket indices
 *
 *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
 */
double FrozenCoreRDMBuilder::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const VectorX<double>& x) const {
    throw std::runtime_error ("FrozenCoreRDMBuilder::calculateElement(std::vector<size_t>, std::vector<size_t>, VectorX<double>): calculateElement is not implemented for FrozenCoreCI RDMs");
};


}  // namespace GQCG
