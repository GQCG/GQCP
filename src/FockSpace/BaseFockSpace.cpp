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
#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/FockSpace.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "FockSpace/SelectedFockSpace.hpp"


namespace GQCP {

/*
 * PROTECTED CONSTRUCTORS
 */

/**
 *  @param K        the number of orbitals
 *  @param dim      the dimension of the Fock space
 */
BaseFockSpace::BaseFockSpace(size_t K, size_t dim) :
    K (K),
    dim (dim)
{}


std::shared_ptr<BaseFockSpace> BaseFockSpace::HeapFockSpace(const BaseFockSpace& fock_space) {

    if (fock_space.is_on_heap){
        //std::shared_ptr<BaseFockSpace> f (&fock_space);
        //return f;
        std::invalid_argument("Passed address of a shared heap fock space");
    }

    std::shared_ptr<BaseFockSpace> fock_space_ptr;

    switch (fock_space.get_type()){

        case FockSpaceType::FockSpace: {
            fock_space_ptr = std::make_shared<FockSpace>(FockSpace(dynamic_cast<const FockSpace&>(fock_space)));
            break;
        }

        case FockSpaceType::ProductFockSpace: {
            fock_space_ptr = std::make_shared<ProductFockSpace>(ProductFockSpace(dynamic_cast<const ProductFockSpace&>(fock_space)));
            break;
        }

        case FockSpaceType::SelectedFockSpace: {
            fock_space_ptr = std::make_shared<SelectedFockSpace>(SelectedFockSpace(dynamic_cast<const SelectedFockSpace&>(fock_space)));
            break;
        }
    }

    fock_space_ptr->is_on_heap = true;
    return fock_space_ptr;
}

/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
BaseFockSpace::~BaseFockSpace() {}



/*
 *  PUBLIC
 */

/**
 *  @return the coefficient vector for the Hartree-Fock wave function (i.e. the 'first' ONV/Slater determinant)
 */
Eigen::VectorXd BaseFockSpace::HartreeFockExpansion() const {
    Eigen::VectorXd expansion = Eigen::VectorXd::Zero(this->dim);
    expansion(0) = 1;  // first configuration is position 0 (conventional ordering of the Fock space)
    return expansion;
}


/**
 *  @return a random normalized coefficient vector, with coefficients uniformly distributed in [-1, 1]
 */
Eigen::VectorXd BaseFockSpace::randomExpansion() const {
    Eigen::VectorXd random = Eigen::VectorXd::Random(this->dim);
    random.normalize();
    return random;
}


/**
 *  @return a constant normalized coefficients vector (i.e. all the coefficients are equal)
 */
Eigen::VectorXd BaseFockSpace::constantExpansion() const {
    Eigen::VectorXd constant = Eigen::VectorXd::Ones(this->dim);
    constant.normalize();
    return constant;
}



}  // namespace GQCP
