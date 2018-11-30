// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#include "FockSpace/ProductFockSpace.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param K            the number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 */
ProductFockSpace::ProductFockSpace(size_t K, size_t N_alpha, size_t N_beta) :
        BaseFockSpace(K, ProductFockSpace::calculateDimension(K, N_alpha, N_beta)),
        fock_space_alpha (FockSpace(K, N_alpha)),
        fock_space_beta (FockSpace(K, N_beta)),
        N_alpha (N_alpha),
        N_beta (N_beta)
{}


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param K            the number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 *
 *  @return the dimension of the product Fock space
 */
size_t ProductFockSpace::calculateDimension(size_t K, size_t N_alpha, size_t N_beta) {
    size_t alpha_dim = FockSpace::calculateDimension(K, N_alpha);
    size_t beta_dim = FockSpace::calculateDimension(K, N_beta);
    return boost::numeric::converter<double, size_t>::convert(beta_dim * alpha_dim);
}


}  // namespace GQCP
