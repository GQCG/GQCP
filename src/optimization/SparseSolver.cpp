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
#include "SparseSolver.hpp"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>



namespace numopt {
namespace eigenproblem {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param dim                      the dimension of the matrix
 *  @param sparse_solver_options    the options to be used for the sparse eigenproblem algorithm
 */
SparseSolver::SparseSolver(size_t dim, const SparseSolverOptions& sparse_solver_options) :
    BaseMatrixSolver (dim, sparse_solver_options.number_of_requested_eigenpairs),
    matrix (Eigen::SparseMatrix<double> (this->dim, this->dim))  // Eigen::Sparse is always initiated to zeros
{}


/**
 *  @param dim                                  the dimension of the matrix
 *  @param number_of_requested_eigenpairs       the number of eigenpairs the eigensolver should find
 */
SparseSolver::SparseSolver(size_t dim, size_t number_of_requested_eigenpairs) :
        BaseMatrixSolver(dim, number_of_requested_eigenpairs),
        matrix (Eigen::SparseMatrix<double> (this->dim, this->dim))  // Eigen::Sparse is always initiated to zeros
{}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Solve the sparse eigenvalue problem
 *
 *  If successful, it sets
 *      - _is_solved to true
 *      - the number of requested eigenpairs
 */
void SparseSolver::solve() {

    // Solve the sparse eigenvalue problem of the Hamiltonian matrix
    Spectra::SparseSymMatProd<double> matrixVectorProduct (this->matrix);

    // Request the number of eigenpairs, and use at least 2 more Ritz pairs than requested eigenvalues)
    Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, Spectra::SparseSymMatProd<double>> spectra_sparse_eigensolver (&matrixVectorProduct, static_cast<int>(this->number_of_requested_eigenpairs), static_cast<int>(this->number_of_requested_eigenpairs) + 2);
    spectra_sparse_eigensolver.init();
    spectra_sparse_eigensolver.compute();

    // Set number of requested eigenpairs

    // Set the eigenvalue and eigenvector as the lowest-energy eigenpair. We can use increasing indices because
    // we have specified Spectra::SMALLEST_ALGE, which selects eigenvalues with smallest algebraic value
    if (spectra_sparse_eigensolver.info() == Spectra::SUCCESSFUL) {
        this->_is_solved = true;

        for (size_t i = 0; i < this->number_of_requested_eigenpairs; i++) {
            double eigenvalue = spectra_sparse_eigensolver.eigenvalues()(i);
            Eigen::VectorXd eigenvector = spectra_sparse_eigensolver.eigenvectors().col(i);

            this->eigenpairs[i] = Eigenpair(eigenvalue, eigenvector);
        }
    } else {  // if Spectra was not successful
        throw std::runtime_error("Spectra could not solve the sparse eigenvalue problem.");
    }
}


/**
 *  @param value        the value to be added
 *  @param index1       the first index of the matrix
 *  @param index2       the second index of the matrix
 *
 *  Add the value to the matrix at (index1, index2)
 */
void SparseSolver::addToMatrix(double value, size_t index1, size_t index2) {
    this->matrix.coeffRef(index1,index2) += value;
}


}  // namespace eigenproblem
}  // namespace numopt
