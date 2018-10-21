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
#ifndef GQCP_BASEFOCKSPACE_HPP
#define GQCP_BASEFOCKSPACE_HPP


#include "ONV.hpp"
#include "FockSpace/FockSpaceType.hpp"
#include "common.hpp"



namespace GQCP {


/**
 *  A base class for the Fock space
 *  Interfacing requires the Fock space to generate an ONV from a given address
 *  transform a given ONV into the next ONV (in the full or selected space)
 *  and retrieve the address of a given ONV in the space
 */
class BaseFockSpace {
protected:
    size_t K;  // number of spatial orbitals
    size_t dim;  // dimension of the Fock space


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor given a @param K and @param dim
     */
    explicit BaseFockSpace(size_t K, size_t dim);
    BaseFockSpace() = default;

public:
    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseFockSpace() = 0;


    // GETTERS
    size_t get_dimension() const { return dim; }
    size_t get_K() const { return K; }
    virtual FockSpaceType get_type() const = 0;

    // PUBLIC METHODS
    /**
     *  Creates a Hartree-Fock coefficient expansion (single Slater expansion of the first configuration in the Fock space)
     */
    Eigen::VectorXd HartreeFockExpansion();
};


}  // namespace GQCP


#endif  // GQCP_BASEFOCKSPACE_HPP
