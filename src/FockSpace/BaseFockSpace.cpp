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
Eigen::VectorXd BaseFockSpace::HartreeFockExpansion() {
    Eigen::VectorXd expansion = Eigen::VectorXd::Zero(this->dim);
    expansion(0) = 1;  // first configuration is position 0 (conventional ordering of the Fock space)
    return expansion;
}


/**
 *  @return a random normalized coefficient vector, with coefficients uniformly distributed in [-1, 1]
 */
Eigen::VectorXd BaseFockSpace::randomExpansion() {
    Eigen::VectorXd random = Eigen::VectorXd::Random(this->dim);
    random.normalize();
    return random;
}


}  // namespace GQCP
