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
#include "DenseSolver.hpp"

#include <iostream>



namespace numopt {
namespace eigenproblem {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param dim                      the dimension of the matrix
 *  @param dense_solver_options     the options to be used for the dense eigenproblem algorithm
 */
DenseSolver::DenseSolver(size_t dim, const DenseSolverOptions& dense_solver_options) :
    BaseMatrixSolver (dim, dense_solver_options.number_of_requested_eigenpairs),
    matrix (Eigen::MatrixXd::Zero(this->dim, this->dim))
{}


/**
 *  @param matrix                   the full dense representation of the matrix
 *  @param dense_solver_options     the options to be used for the dense eigenproblem algorithm
 */
DenseSolver::DenseSolver(const Eigen::MatrixXd& matrix, const DenseSolverOptions& dense_solver_options) :
        BaseMatrixSolver (matrix.cols(), dense_solver_options.number_of_requested_eigenpairs),
        matrix (matrix)
{}

    
/**
 *  @param dim                                  the dimension of the matrix
 *  @param number_of_requested_eigenpairs       the number of eigenpairs the eigensolver should find
 */
DenseSolver::DenseSolver(size_t dim, size_t number_of_requested_eigenpairs) :
        BaseMatrixSolver (dim, number_of_requested_eigenpairs),
        matrix (Eigen::MatrixXd::Zero(this->dim, this->dim))
{}


/**
 *  @param matrix                   the full dense representation of the matrix
 *  @param number_of_requested_eigenpairs       the number of eigenpairs the eigensolver should find
 */
DenseSolver::DenseSolver(Eigen::MatrixXd matrix, size_t number_of_requested_eigenpairs) :
        BaseMatrixSolver (matrix.cols(), number_of_requested_eigenpairs),
        matrix (matrix)
{}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Solve the dense eigenvalue problem
 *
 *  If successful, it sets
 *      - @member is_solved to true
 *      - the number of requested eigenpairs in @member eigenpairs
 */
void DenseSolver::solve() {

    // Solve the dense eigenvalue problem of the Hamiltonian matrix.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> self_adjoint_eigensolver (this->matrix);  // gives the eigenvalues (and corresponding eigenvectors) in ascending order


    // Set number of requested eigenpairs
    this->_is_solved = true;
    for (size_t i = 0; i < this->number_of_requested_eigenpairs; i++) {
        double eigenvalue = self_adjoint_eigensolver.eigenvalues()(i);
        Eigen::VectorXd eigenvector = self_adjoint_eigensolver.eigenvectors().col(i);

        this->eigenpairs[i] = Eigenpair(eigenvalue, eigenvector);
    }
}


/**
 *  @param value        the value to be added
 *  @param index1       the first index of the matrix
 *  @param index2       the second index of the matrix
 *
 *  Add the value to the matrix at (index1, index2)
 */
void DenseSolver::addToMatrix(double value, size_t index1, size_t index2) {
    this->matrix(index1, index2) += value;
}


}  // namespace eigenproblem
}  // namespace numopt
