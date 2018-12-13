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
#ifndef NUMOPT_BASEMATRIXSOLVER_HPP
#define NUMOPT_BASEMATRIXSOLVER_HPP



#include "BaseEigenproblemSolver.hpp"



namespace numopt {
namespace eigenproblem {


/**
 *  A base class for eigenproblem solvers in which the whole matrix is supplied
 */
class BaseMatrixSolver : public numopt::eigenproblem::BaseEigenproblemSolver {
public:
    // CONSTRUCTOR
    /**
     *  @param dim                                  the dimension of the eigenvalue problem
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
     */
    explicit BaseMatrixSolver(size_t dim, size_t number_of_requested_eigenpairs = 1);


    // DESTRUCTOR
    ~BaseMatrixSolver() override = default;


    // PUBLIC PURE VIRTUAL METHODS
    /**
     *  @param value        the value to be added
     *  @param index1       the first index of the matrix
     *  @param index2       the second index of the matrix
     *
     *  Add the value to the matrix at (index1, index2)
     */
    virtual void addToMatrix(double value, size_t index1, size_t index2) = 0;
};


}  // namespace eigenproblem
}  // namespace numopt



#endif  // NUMOPT_BASEMATRIXSOLVER_HPP
