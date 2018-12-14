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
#include "optimization/BaseEigenproblemSolver.hpp"
#include <iostream>


namespace GQCP {


/*
 *  PROTECTED CONSTRUCTORS
 */

/**
 *  @param dim                                  the dimension of the vector space associated to the eigenvalue problem
 *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
 */
BaseEigenproblemSolver::BaseEigenproblemSolver(size_t dim, size_t number_of_requested_eigenpairs) :
    dim (dim),
    number_of_requested_eigenpairs (number_of_requested_eigenpairs)
{
    eigenpairs.reserve(this->number_of_requested_eigenpairs);
}



/*
 *  GETTERS
 */

std::vector<const Eigenpair> BaseEigenproblemSolver::get_eigenpairs() const {

    if (this->_is_solved) {
        return this->eigenpairs;
    } else {
        throw std::logic_error("The eigenvalue problem hasn't been solved yet and you are trying to get the eigenpairs.");
    }
}


/**
 *  @param i        the index for the i-th lowest eigenpair
 *
 *  @return the i-th lowest eigenpair
 */
const Eigenpair& BaseEigenproblemSolver::get_eigenpair(size_t i) const {

    return this->get_eigenpairs()[i];
}


/**
 *  @param i        the index for the i-th lowest eigenvalue
 *
 *  @return the i-th lowest eigenvalue
 */
double BaseEigenproblemSolver::get_eigenvalue(size_t i) const {

    return this->get_eigenpair(i).get_eigenvalue();
}


/**
 *  @param i        the index for the i-th lowest eigenvector
 *
 *  @return the i-th lowest eigenvector
 */
const Eigen::VectorXd& BaseEigenproblemSolver::get_eigenvector(size_t i) const {

    return this->get_eigenpair(i).get_eigenvector();
}



}  // namespace GQCP
