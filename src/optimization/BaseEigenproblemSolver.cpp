// This file is part of GQCG-numopt.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-numopt is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-numopt is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-numopt.  If not, see <http://www.gnu.org/licenses/>.
#include "BaseEigenproblemSolver.hpp"
#include <iostream>


namespace numopt {
namespace eigenproblem {


/*
 *  PROTECTED CONSTRUCTORS
 */

/**
 *  @param dim                                  the dimension of the vector space associated to the eigenvalue problem
 *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
 */
BaseEigenproblemSolver::BaseEigenproblemSolver(size_t dim, size_t number_of_requested_eigenpairs) :
    dim (dim),
    number_of_requested_eigenpairs (number_of_requested_eigenpairs),
    eigenpairs (std::vector<numopt::eigenproblem::Eigenpair> (this->number_of_requested_eigenpairs, numopt::eigenproblem::Eigenpair(this->dim)))
{}



/*
 *  GETTERS
 */
std::vector<numopt::eigenproblem::Eigenpair> BaseEigenproblemSolver::get_eigenpairs() const {

    if (this->_is_solved) {
        return this->eigenpairs;
    } else {
        throw std::logic_error("The eigenvalue problem hasn't been solved yet and you are trying to get the eigenpairs.");
    }
}


numopt::eigenproblem::Eigenpair BaseEigenproblemSolver::get_lowest_eigenpair() const {

    if (this->_is_solved) {
        return this->eigenpairs[0];  // the eigenpairs are sorted by increasing eigenvalue
    } else {
        throw std::logic_error("The eigenvalue problem hasn't been solved yet and you are trying to get the lowest eigenpair.");
    }
}


/**
 *  @param i        the index for the i-th lowest eigenpair
 *
 *  @return the i-th lowest eigenpair
 */
numopt::eigenproblem::Eigenpair BaseEigenproblemSolver::get_eigenpair(size_t i) const {

    if (this->_is_solved) {
        return this->eigenpairs[i];  // the eigenpairs are sorted by increasing eigenvalue
    } else {
        throw std::logic_error("The eigenvalue problem hasn't been solved yet and you are trying to get the i-th lowest eigenpair.");
    }
}


double BaseEigenproblemSolver::get_lowest_eigenvalue() const {

    if (this->_is_solved) {
        return this->eigenpairs[0].get_eigenvalue();
    } else {
        throw std::logic_error("The eigenvalue problem hasn't been solved yet and you are trying to get the lowest eigenvalue.");
    }
}


/**
 *  Special shortcut getter for the lowest eigenvalue: will be deprecated in the next major release
 */
double BaseEigenproblemSolver::get_eigenvalue() const {

    return this->get_lowest_eigenvalue();
}


Eigen::VectorXd BaseEigenproblemSolver::get_lowest_eigenvector() const {

    if (this->_is_solved) {
        return this->eigenpairs[0].get_eigenvector();
    } else {
        throw std::logic_error("The eigenvalue problem hasn't been solved yet and you are trying to get the lowest eigenvector.");
    }
}


/**
 *  @param index        the index of an element of the lowest eigenvector
 *
 *  @return the element at the given index of the lowest eigenvector
 */
double BaseEigenproblemSolver::get_lowest_eigenvector(size_t index) const {

    return this->get_lowest_eigenvector()[index];
}


/**
 *  Special shortcut getter for the eigenvector corresponding to the lowest eigenvalue: will be deprecated in the next major release
 */
Eigen::VectorXd BaseEigenproblemSolver::get_eigenvector() const {

    return this->get_lowest_eigenvector();
}


/**
 *  @param index        the index of an element of the lowest eigenvector
 *
 *  @return the value at the index of the eigenvector corresponding to the lowest eigenvalue: will be deprecated in the next major release
 */
double BaseEigenproblemSolver::get_eigenvector(size_t index) const {

    return this->get_lowest_eigenvector(index);
}



}  // namespace eigenproblem
}  // namespace numopt
