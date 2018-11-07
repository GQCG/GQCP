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
 *  A base class for the representation of a Fock space
 */
class BaseFockSpace {
protected:
    size_t K;  // number of spatial orbitals
    size_t dim;  // dimension of the Fock space


    // PROTECTED CONSTRUCTORS
    BaseFockSpace() = default;
    /**
     *  @param K        the number of orbitals
     *  @param dim      the dimension of the Fock space
     */
    explicit BaseFockSpace(size_t K, size_t dim);

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
     *  @return the coefficient vector for the Hartree-Fock wave function (i.e. the 'first' ONV/Slater determinant)
     */
    Eigen::VectorXd HartreeFockExpansion();

    /**
     *  @return a random normalized coefficient vector, with coefficients uniformly distributed in [-1, 1]
     */
    Eigen::VectorXd randomExpansion();
};


}  // namespace GQCP


#endif  // GQCP_BASEFOCKSPACE_HPP
