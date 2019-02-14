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
#ifndef GQCP_SPARSESOLVER_HPP
#define GQCP_SPARSESOLVER_HPP



#include <Eigen/Sparse>

#include "optimization/BaseMatrixSolver.hpp"
#include "optimization/EigenproblemSolverOptions.hpp"



namespace GQCP {


/**
 *  An eigenproblem solver that stores a sparse representation of the matrix
 */
class SparseSolver : public BaseMatrixSolver {
private:
    Eigen::SparseMatrix<double> matrix;


public:
    // CONSTRUCTORS
    /**
     *  Constructor that sets the sparse matrix to zeros
     *
     *  @param dim                                  the dimension of the matrix
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the eigensolver should find
     */
    SparseSolver(size_t dim, size_t number_of_requested_eigenpairs = 1);

    /**
     *  Constructor that sets the sparse matrix to zeros
     *
     *  @param dim                      the dimension of the matrix
     *  @param sparse_solver_options    the options to be used for the sparse eigenproblem algorithm
     */
    SparseSolver(size_t dim, const SparseSolverOptions& sparse_solver_options);


    // DESTRUCTOR
    ~SparseSolver() override = default;


    // GETTERS
    const Eigen::SparseMatrix<double>& get_matrix() const { return this->matrix; }


    // PUBLIC OVERRIDDEN METHODS
    /**
     *  Solve the sparse eigenvalue problem
     *
     *  If successful, it sets
     *      - _is_solved to true
     *      - the number of requested eigenpairs
     */
    void solve() override;

    /**
     *  @param value        the value to be added
     *  @param index1       the first index of the matrix
     *  @param index2       the second index of the matrix
     *
     *  Add the value to the matrix at (index1, index2)
     */
    void addToMatrix(double value, size_t index1, size_t index2) override;
};


}  // namespace GQCP



#endif  // GQCP_SPARSESOLVER_HPP
