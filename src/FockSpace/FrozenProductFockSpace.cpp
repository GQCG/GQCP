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
#include "FockSpace/FrozenProductFockSpace.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param K            the number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 *  @param X            the number of frozen orbitals
 */
FrozenProductFockSpace::FrozenProductFockSpace(size_t K, size_t N_alpha, size_t N_beta, size_t X) :
        BaseFockSpace(K, ProductFockSpace::calculateDimension(K-X, N_alpha-X, N_beta-X)),
        fock_space_alpha (FrozenFockSpace(K, N_alpha, X)),
        fock_space_beta (FrozenFockSpace(K, N_beta, X)),
        fock_space (K-X, N_alpha-X, N_beta-X),
        X (X)
{}


/**
 *  @param fock_space       non-frozen sub product Fock space
 *  @param X                the number of frozen orbitals
 */
FrozenProductFockSpace::FrozenProductFockSpace(const ProductFockSpace& fock_space, size_t X) :
    FrozenProductFockSpace(fock_space.get_K(), fock_space.get_N_alpha(), fock_space.get_N_beta(), X)
{}

}  // namespace GQCP
