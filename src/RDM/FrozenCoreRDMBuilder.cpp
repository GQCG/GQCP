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
#include "RDM/BaseRDMBuilder.hpp"


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
OneRDMs calculate1RDMs(const Eigen::VectorXd& x) const {

    auto K = this->rdm_builder->get_fock_space().get_K();

    Eigen::MatrixXd D_aa = Eigen::MatrixXd::Zero(K, K);
    Eigen::MatrixXd D_bb = Eigen::MatrixXd::Zero(K, K);

    auto Kn = K - X;

    for (size_t i = 0; i < X; i++) {
        D_aa(i,i) = 1;
        D_bb(i,i) = 1;
    }

    OneRDMs sub_1rdms = this->rdm_builder->calculate1RDMs(x);

    D_aa.block(X, X, Kn, Kn) += sub_1rdms.one_rdm_aa.get_matrix_representation();
    D_bb.block(X, X, Kn, Kn) += sub_1rdms.one_rdm_bb.get_matrix_representation();

    OneRDM one_rdm_aa (D_aa);
    OneRDM one_rdm_bb (D_bb);

    return OneRDMs (one_rdm_aa, one_rdm_bb);
};

/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 2-RDMs given a coefficient vector
 */
TwoRDMs calculate2RDMs(const Eigen::VectorXd& x) const {

    auto K = this->rdm_builder->get_fock_space().get_K();
    OneRDMs one_rdms = this->calculate1RDMs(x);
    TwoRDMs sub_2rdms = this->rdm_builder->calculate2RDMs(x);


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
double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const Eigen::VectorXd& x) const {
    throw std::runtime_error ("calculateElement is not implemented for FrozenCore RDMs");
};

}  // namespace GQCG
