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
#ifndef NUMOPT_BASEEIGENVALUESOLVER_HPP
#define NUMOPT_BASEEIGENVALUESOLVER_HPP


#include "Eigenpair.hpp"

#include <cstddef>
#include <Eigen/Dense>
#include <vector>


namespace numopt {
namespace eigenproblem {


/**
 *  A base class for the implementation of eigenvalue problem solvers.
 *
 *  Derived classes should implement
 *      - get_diagonal()
 */
class BaseEigenproblemSolver {
protected:
    const size_t dim;  // the dimension of the vector space associated to the eigenvalue problem
    const size_t number_of_requested_eigenpairs;

    bool _is_solved = false;
    std::vector<numopt::eigenproblem::Eigenpair> eigenpairs;  // a collection of the eigenpairs of the eigenproblem
                                                              // the eigenpairs are sorted with increasing eigenvalue


    // PROTECTED CONSTRUCTORS
    /**
     *  @param dim                                  the dimension of the vector space associated to the eigenvalue problem
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
     */
    explicit BaseEigenproblemSolver(size_t dim, size_t number_of_requested_eigenpairs = 1);


public:
    // DESTRUCTOR
    virtual ~BaseEigenproblemSolver() = default;


    // GETTERS
    /**
     *  @return the diagonal of the matrix representation of the matrix whose eigenproblem is being solved
     */
    virtual Eigen::VectorXd get_diagonal() = 0;
    bool is_solved() const { return this->_is_solved; }


    // GETTERS - EIGENPAIR
    std::vector<numopt::eigenproblem::Eigenpair> get_eigenpairs() const;
    numopt::eigenproblem::Eigenpair get_lowest_eigenpair() const;
    /**
     *  @param i        the index for the i-th lowest eigenpair
     *
     *  @return the i-th lowest eigenpair
     */
    numopt::eigenproblem::Eigenpair get_eigenpair(size_t i) const;


    // GETTERS - EIGENVALUE
    double get_lowest_eigenvalue() const;
    /**
     *  Special shortcut getter for the lowest eigenvalue: will be deprecated in the next major release
     */
    double get_eigenvalue() const;


    // GETTERS - EIGENVECTOR
    Eigen::VectorXd get_lowest_eigenvector() const;
    /**
     *  @param index        the index of an element of the lowest eigenvector
     *
     *  @return the element at the given index of the lowest eigenvector
     */
    double get_lowest_eigenvector(size_t index) const;
    /**
     *  Special shortcut getter for the eigenvector corresponding to the lowest eigenvalue: will be deprecated in the next major release
     */
    Eigen::VectorXd get_eigenvector() const;
    /**
     *  @param index        the index of an element of the lowest eigenvector
     *
     *  @return the value at the index of the eigenvector corresponding to the lowest eigenvalue: will be deprecated in the next major release
     */
    double get_eigenvector(size_t index) const;


    // PUBLIC PURE VIRTUAL METHODS
    /**
     *  Solve the eigenvalue problem associated to the eigenproblem solver
     *
     *  If successful, it sets
     *      - @member is_solved to true
     *      - the number of requested eigenpairs in @member eigenpairs
     */
    virtual void solve() = 0;
};


}  // namespace eigenproblem
}  // namespace numopt





#endif  // NUMOPT_BASEEIGENVALUESOLVER_HPP
