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
#ifndef GQCP_FROZENPRODUCTFOCKSPACE_HPP
#define GQCP_FROZENPRODUCTFOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "FockSpace/FrozenFockSpace.hpp"


namespace GQCP {


/**
 *  A class that represents the product of two frozen Fock spaces (alpha and beta).
 */
class FrozenProductFockSpace: public BaseFockSpace {
private:
    size_t X;  // number of frozen orbitals/electrons

    FrozenFockSpace frozen_fock_space_alpha;
    FrozenFockSpace frozen_fock_space_beta;

    ProductFockSpace active_product_fock_space;  // active (non-frozen) product Fock space containing only the active electrons (N_alpha-X, N_beta-X) and orbitals (K-X)

public:
    // CONSTRUCTORS
    /**
     *  @param K            the total number of orbitals (equal for alpha and beta)
     *  @param N_alpha      the total number of alpha electrons
     *  @param N_beta       the total number of beta electrons
     *  @param X            the number of frozen orbitals and electrons (equal for alpha and beta)
     */
    FrozenProductFockSpace(size_t K, size_t N_alpha, size_t N_beta, size_t X);

    /**
     *  @param fock_space       (to be frozen) full product Fock space
     *  @param X                the number of frozen orbitals and electrons (equal for alpha and beta)
     */
    FrozenProductFockSpace(const ProductFockSpace& fock_space, size_t X);

    // GETTERS
    size_t get_N_alpha() const { return this->frozen_fock_space_alpha.get_N(); }
    size_t get_N_beta() const { return this->frozen_fock_space_beta.get_N(); }
    size_t get_number_of_frozen_orbitals() const { return this->X;}

    const FrozenFockSpace& get_frozen_fock_space_alpha() const { return this->frozen_fock_space_alpha; }
    const FrozenFockSpace& get_frozen_fock_space_beta() const { return this->frozen_fock_space_beta; }

    const ProductFockSpace& get_active_product_fock_space() const { return this->active_product_fock_space; }
    FockSpaceType get_type() const override { return FockSpaceType::FrozenProductFockSpace; }
};


}  // namespace GQCP


#endif  // GQCP_FOCKSPACE_HPP
