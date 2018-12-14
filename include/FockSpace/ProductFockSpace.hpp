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
#ifndef GQCP_PRODUCTFOCKSPACE_HPP
#define GQCP_PRODUCTFOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/FockSpace.hpp"


namespace GQCP {


/**
 *  A class that represents the product of two full Fock spaces (alpha and beta).
 */
class ProductFockSpace: public BaseFockSpace {
private:
    size_t N_alpha;  // number of alpha electrons
    size_t N_beta;  // number of beta electrons

    FockSpace fock_space_alpha;
    FockSpace fock_space_beta;


public:
    // CONSTRUCTORS
    /**
     *  @param K            the number of orbitals (equal for alpha and beta)
     *  @param N_alpha      the number of alpha electrons
     *  @param N_beta       the number of beta electrons
     */
    ProductFockSpace(size_t K, size_t N_alpha, size_t N_beta);


    // DESTRUCTORS
    ~ProductFockSpace() override = default;


    // GETTERS
    size_t get_N_alpha() const { return this->N_alpha; }
    size_t get_N_beta() const { return this->N_beta; }
    const FockSpace& get_fock_space_alpha() const { return this->fock_space_alpha; }
    const FockSpace& get_fock_space_beta() const { return this->fock_space_beta; }
    FockSpaceType get_type() const override { return FockSpaceType::ProductFockSpace; }


    // STATIC PUBLIC METHODS
    /**
     *  @param K            the number of orbitals (equal for alpha and beta)
     *  @param N_alpha      the number of alpha electrons
     *  @param N_beta       the number of beta electrons
     *
     *  @return the dimension of the product Fock space
     */
    static size_t calculateDimension(size_t K, size_t N_alpha, size_t N_beta);
};


}  // namespace GQCP


#endif  // GQCP_FOCKSPACE_HPP
