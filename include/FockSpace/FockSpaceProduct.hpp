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
#ifndef GQCP_FOCKSPACEPRODUCT_HPP
#define GQCP_FOCKSPACEPRODUCT_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/FockSpace.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCP {


/**
 *  The product of two Fock spaces for a given set of orbitals and number of alpha and beta electrons.
 */
class FockSpaceProduct: public GQCP::BaseFockSpace {
private:
    const size_t N_alpha;  // number of alpha electrons
    const size_t N_beta;  // number of beta electrons

    FockSpace fock_space_alpha;
    FockSpace fock_space_beta;


public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param K (spatial orbitals), N_alpha and N_beta (electrons)
     *  on which the dimensions of the Fock space are based
     */
    FockSpaceProduct(size_t K, size_t N_alpha, size_t N_beta);


    // DESTRUCTORS
    ~FockSpaceProduct() override = default;


    // GETTERS
    size_t get_N_alpha() const { return this->N_alpha; }
    size_t get_N_beta() const { return this->N_beta; }
    FockSpace get_fock_space_alpha() const { return this->fock_space_alpha; }
    FockSpace get_fock_space_beta() const { return this->fock_space_beta; }
    FockSpaceType get_fock_space_type() const override { return FockSpaceType::FockSpaceProduct; }


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K
     *  and a number of alpha and beta electrons @param N_alpha, N_beta,
     *  @return the dimension of the Fock space
     */
    static size_t calculateDimension(size_t K, size_t N_alpha, size_t N_beta);
};


}  // namespace GQCP


#endif  // GQCP_FOCKSPACE_HPP
