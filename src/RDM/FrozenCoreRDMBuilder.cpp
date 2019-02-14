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
FrozenCoreRDMBuilder::FrozenCoreRDMBuilder(std::shared_ptr<BaseRDMBuilder> rdm_builder, size_t X) :
        BaseRDMBuilder(),
        rdm_builder (std::move(rdm_builder)),
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
OneRDMs FrozenCoreRDMBuilder::calculate1RDMs(const Eigen::VectorXd& x) const {

    auto K = this->rdm_builder->get_fock_space()->get_K() + this->X;

    Eigen::MatrixXd D_aa = Eigen::MatrixXd::Zero(K, K);
    Eigen::MatrixXd D_bb = Eigen::MatrixXd::Zero(K, K);

    auto Kn = K - this->X;

    for (size_t i = 0; i < this->X; i++) {
        D_aa(i,i) = 1;
        D_bb(i,i) = 1;
    }

    // Retrieve 1RDMs from non-frozen sub-space
    OneRDMs sub_1rdms = this->rdm_builder->calculate1RDMs(x);

    // Incorporate sub rdms into total rdms
    D_aa.block(this->X, this->X, Kn, Kn) += sub_1rdms.one_rdm_aa.get_matrix_representation();
    D_bb.block(this->X, this->X, Kn, Kn) += sub_1rdms.one_rdm_bb.get_matrix_representation();

    OneRDM one_rdm_aa (D_aa);
    OneRDM one_rdm_bb (D_bb);

    return OneRDMs (one_rdm_aa, one_rdm_bb);
};

/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 2-RDMs given a coefficient vector
 */
TwoRDMs FrozenCoreRDMBuilder::calculate2RDMs(const Eigen::VectorXd& x) const {

    auto K = this->rdm_builder->get_fock_space()->get_K() + this->X;

    Eigen::Tensor<double, 4> d_aaaa (K,K,K,K);
    d_aaaa.setZero();
    Eigen::Tensor<double, 4> d_aabb (K,K,K,K);
    d_aabb.setZero();
    Eigen::Tensor<double, 4> d_bbaa (K,K,K,K);
    d_bbaa.setZero();
    Eigen::Tensor<double, 4> d_bbbb (K,K,K,K);
    d_bbbb.setZero();

    // Retrieve sub one RDMs to implement the frozen RDM formulas
    OneRDMs one_rdms = this->rdm_builder->calculate1RDMs(x);
    auto D_aa = one_rdms.one_rdm_aa.get_matrix_representation();
    auto D_bb = one_rdms.one_rdm_bb.get_matrix_representation();

    // Implement frozen RDM formulas
    for (size_t p = 0; p < this->X; p++) {

        GQCP::tensorBlockAddition<0,1>(d_aaaa, D_aa, this->X, this->X, p, p);
        GQCP::tensorBlockAddition<2,3>(d_aaaa, D_aa, p, p, this->X, this->X);
        GQCP::tensorBlockAddition<2,1>(d_aaaa, -1*D_aa, p, this->X, this->X, p);
        GQCP::tensorBlockAddition<3,0>(d_aaaa, -1*D_aa, this->X, p, p, this->X);

        GQCP::tensorBlockAddition<0,1>(d_bbbb, D_bb, this->X, this->X, p, p);
        GQCP::tensorBlockAddition<2,3>(d_bbbb, D_bb, p, p, this->X, this->X);
        GQCP::tensorBlockAddition<2,1>(d_bbbb, -1*D_bb, p, this->X, this->X, p);
        GQCP::tensorBlockAddition<3,0>(d_bbbb, -1*D_bb, this->X, p, p, this->X);

        GQCP::tensorBlockAddition<2,3>(d_aabb, D_bb, p, p, this->X, this->X);
        GQCP::tensorBlockAddition<0,1>(d_aabb, D_aa, this->X, this->X, p, p);

        GQCP::tensorBlockAddition<2,3>(d_bbaa, D_aa, p, p, this->X, this->X);
        GQCP::tensorBlockAddition<0,1>(d_bbaa, D_bb, this->X, this->X, p, p);

        d_bbaa(p,p,p,p) = 1;
        d_aabb(p,p,p,p) = 1;

        for (size_t q = p +1; q < this->X; q++) {

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

    // retrieve non frozen sub-space 2RDMs
    TwoRDMs sub_2rdms = this->rdm_builder->calculate2RDMs(x);

    // incorporate into the total 2RDMs
    GQCP::tensorBlockAddition(d_aaaa, sub_2rdms.two_rdm_aaaa.get_matrix_representation(), this->X, this->X, this->X, this->X);
    GQCP::tensorBlockAddition(d_bbbb, sub_2rdms.two_rdm_bbbb.get_matrix_representation(), this->X, this->X, this->X, this->X);
    GQCP::tensorBlockAddition(d_aabb, sub_2rdms.two_rdm_aabb.get_matrix_representation(), this->X, this->X, this->X, this->X);
    GQCP::tensorBlockAddition(d_bbaa, sub_2rdms.two_rdm_bbaa.get_matrix_representation(), this->X, this->X, this->X, this->X);

    TwoRDM two_rdm_aaaa (d_aaaa);
    TwoRDM two_rdm_aabb (d_aabb);
    TwoRDM two_rdm_bbaa (d_bbaa);
    TwoRDM two_rdm_bbbb (d_bbbb);

    return TwoRDMs (two_rdm_aaaa, two_rdm_aabb, two_rdm_bbaa, two_rdm_bbbb);
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
double FrozenCoreRDMBuilder::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const Eigen::VectorXd& x) const {
    throw std::runtime_error ("calculateElement is not implemented for FrozenCore RDMs");
};

}  // namespace GQCG
