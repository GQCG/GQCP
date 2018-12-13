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
#include "BaseMatrixSolver.hpp"



namespace numopt {
namespace eigenproblem {



// CONSTRUCTOR
/**
 *  @param dim                                  the dimension of the eigenvalue problem
 *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
 */
BaseMatrixSolver::BaseMatrixSolver(size_t dim, size_t number_of_requested_eigenpairs) :
    BaseEigenproblemSolver(dim, number_of_requested_eigenpairs)
{}



}  // namespace eigenproblem
}  // namespace numopt
