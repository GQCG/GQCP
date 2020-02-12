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
#include "ONVBasis/BaseONVBasis.hpp"

#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"


namespace GQCP {


/*
 * PROTECTED CONSTRUCTORS
 */

/**
 *  @param K        the number of orbitals
 *  @param dim      the dimension of the ONV basis
 */
BaseONVBasis::BaseONVBasis(size_t K, size_t dim) :
    K (K),
    dim (dim)
{}



/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Clones a derived BaseONVBasis instance to the heap memory
 *
 *  @param fock_space     reference to a derived BaseONVBasis instance to be cloned.
 *
 *  @return a shared pointer owning the heap-cloned ONV basis
 */
std::shared_ptr<BaseONVBasis> BaseONVBasis::CloneToHeap(const BaseONVBasis& fock_space) {

    std::shared_ptr<BaseONVBasis> fock_space_ptr;

    switch (fock_space.get_type()) {

        case ONVBasisType::SpinUnresolvedONVBasis: {
            fock_space_ptr = std::make_shared<SpinUnresolvedONVBasis>(SpinUnresolvedONVBasis(dynamic_cast<const SpinUnresolvedONVBasis&>(fock_space)));
            break;
        }

        case ONVBasisType::SpinResolvedONVBasis: {
            fock_space_ptr = std::make_shared<SpinResolvedONVBasis>(SpinResolvedONVBasis(dynamic_cast<const SpinResolvedONVBasis&>(fock_space)));
            break;
        }

        case ONVBasisType::SpinResolvedSelectedONVBasis: {
            fock_space_ptr = std::make_shared<SpinResolvedSelectedONVBasis>(SpinResolvedSelectedONVBasis(dynamic_cast<const SpinResolvedSelectedONVBasis&>(fock_space)));
            break;
        }

        case ONVBasisType::SpinUnresolvedFrozenONVBasis: {
            fock_space_ptr = std::make_shared<SpinUnresolvedFrozenONVBasis>(SpinUnresolvedFrozenONVBasis(dynamic_cast<const SpinUnresolvedFrozenONVBasis&>(fock_space)));
            break;
        }

        case ONVBasisType::SpinResolvedFrozenONVBasis: {
            fock_space_ptr = std::make_shared<SpinResolvedFrozenONVBasis>(SpinResolvedFrozenONVBasis(dynamic_cast<const SpinResolvedFrozenONVBasis&>(fock_space)));
            break;
        }
    }

    return fock_space_ptr;
}



/*
 *  PUBLIC
 */

/**
 *  @return the coefficient vector for the Hartree-Fock wave function (i.e. the 'first' ONV/Slater determinant)
 */
VectorX<double> BaseONVBasis::hartreeFockExpansion() const {
    VectorX<double> expansion = VectorX<double>::Zero(this->dim);
    expansion(0) = 1;  // first configuration is position 0 (conventional ordering of the ONV basis)
    return expansion;
}


/**
 *  @return a random normalized coefficient vector, with coefficients uniformly distributed in [-1, 1]
 */
VectorX<double> BaseONVBasis::randomExpansion() const {
    VectorX<double> random = VectorX<double>::Random(this->dim);
    random.normalize();
    return random;
}


/**
 *  @return a constant normalized coefficients vector (i.e. all the coefficients are equal)
 */
VectorX<double> BaseONVBasis::constantExpansion() const {
    VectorX<double> constant = VectorX<double>::Ones(this->dim);
    constant.normalize();
    return constant;
}


}  // namespace GQCP
