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

#include "ONVBasis/ONVBasis.hpp"
#include "ONVBasis/ProductONVBasis.hpp"
#include "ONVBasis/SelectedONVBasis.hpp"


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

        case ONVBasisType::ONVBasis: {
            fock_space_ptr = std::make_shared<ONVBasis>(ONVBasis(dynamic_cast<const ONVBasis&>(fock_space)));
            break;
        }

        case ONVBasisType::ProductONVBasis: {
            fock_space_ptr = std::make_shared<ProductONVBasis>(ProductONVBasis(dynamic_cast<const ProductONVBasis&>(fock_space)));
            break;
        }

        case ONVBasisType::SelectedONVBasis: {
            fock_space_ptr = std::make_shared<SelectedONVBasis>(SelectedONVBasis(dynamic_cast<const SelectedONVBasis&>(fock_space)));
            break;
        }

        case ONVBasisType::FrozenONVBasis: {
            fock_space_ptr = std::make_shared<FrozenONVBasis>(FrozenONVBasis(dynamic_cast<const FrozenONVBasis&>(fock_space)));
            break;
        }

        case ONVBasisType::FrozenProductONVBasis: {
            fock_space_ptr = std::make_shared<FrozenProductONVBasis>(FrozenProductONVBasis(dynamic_cast<const FrozenProductONVBasis&>(fock_space)));
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
VectorX<double> BaseONVBasis::HartreeFockExpansion() const {
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
