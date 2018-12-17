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
#ifndef GQCP_BASEEIGENPROBLEMSOLVER_HPP
#define GQCP_BASEEIGENPROBLEMSOLVER_HPP


#include "optimization/Eigenpair.hpp"

#include <cstddef>
#include <Eigen/Dense>
#include <vector>


namespace GQCP {


/**
 *  A base class for the implementation of eigenvalue problem solvers
 */
class BaseEigenproblemSolver {
protected:
    size_t dim;  // the dimension of the vector space associated to the eigenvalue problem
    size_t number_of_requested_eigenpairs;

    bool _is_solved = false;
    std::vector<Eigenpair> eigenpairs;  // a collection of the eigenpairs of the eigenproblem, sorted with increasing eigenvalue


    // PROTECTED CONSTRUCTORS
    /**
     *  @param dim                                  the dimension of the vector space associated to the eigenvalue problem
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
     */
    BaseEigenproblemSolver(size_t dim, size_t number_of_requested_eigenpairs = 1);


public:
    // DESTRUCTOR
    virtual ~BaseEigenproblemSolver() = default;


    // GETTERS
    bool is_solved() const { return this->_is_solved; }

    const std::vector<Eigenpair>& get_eigenpairs() const;

    /**
     *  @param i        the index for the i-th lowest eigenpair
     *
     *  @return the i-th lowest eigenpair
     */
    const Eigenpair& get_eigenpair(size_t i = 0) const;

    /**
     *  @param i        the index for the i-th lowest eigenvalue
     *
     *  @return the i-th lowest eigenvalue
     */
    double get_eigenvalue(size_t i = 0) const;

    /**
     *  @param i        the index for the i-th lowest eigenvector
     *
     *  @return the i-th lowest eigenvector
     */
    const Eigen::VectorXd& get_eigenvector(size_t i = 0) const;


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


}  // namespace GQCP





#endif  // GQCP_BASEEIGENPROBLEMSOLVER_HPP
