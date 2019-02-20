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
#ifndef GQCP_FROZENCORERDMBUILDER_HPP
#define GQCP_FROZENCORERDMBUILDER_HPP


#include "RDM/BaseRDMBuilder.hpp"
#include "RDM/RDMs.hpp"


namespace GQCP {


/**
 *  A class capable of calculating 1- and 2-RDMs from wave functions expanded in the full CI product Fock space
 */
class FrozenCoreRDMBuilder : public BaseRDMBuilder {
    size_t X;  // number of frozen orbitals/electrons
    std::shared_ptr<BaseRDMBuilder> rdm_builder;  // active (non-frozen core) RDM builder performing the BaseRDMBuilder interface in the active space with the frozen core CI wave function

public:
    // CONSTRUCTORS
    /**
     *  @param rdm_builder                  shared pointer to active (non-frozen core) RDM builder
     *  @param X                            the number of frozen orbitals
     */
    FrozenCoreRDMBuilder(std::shared_ptr<BaseRDMBuilder> rdm_builder, size_t X);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 1-RDMs given a coefficient vector
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd& x) const override;

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 2-RDMs given a coefficient vector
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd& x) const override;

    /**
     *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
     *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
     *  @param x                the coefficient vector representing the wave function
     *
     *  @return an element of the spin-summed (total) N-RDM, as specified by the given bra and ket indices
     *
     *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
     */
    double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const Eigen::VectorXd& x) const override;
};


}  // namespace GQCP


#endif  // GQCP_FROZENCORERDMBUILDER_HPP
