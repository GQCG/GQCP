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
#pragma once


#include "Mathematical/Optimization/Eigenproblem/Eigenpair.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  An environment used to solve eigenvalue problems for self-adjoint matrices.
 */
class EigenproblemEnvironment {
public:
    VectorFunction<double> matrix_vector_product_function;  // a vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix)
    SquareMatrix<double> A;  // the self-adjoint matrix whose eigenvalue problem should be solved
    VectorX<double> diagonal;  // the diagonal of the matrix
    size_t dimension;  // the dimension of the diagonalization problem (the dimension of one eigenvector)

    VectorX<double> eigenvalues;  // the eigenvalues of the matrix A
    MatrixX<double> eigenvectors;  // the eigenvectors of the matrix A


    SquareMatrix<double> S;  // the "subspace matrix": the projection of the matrix A onto the subspace spanned by the vectors in V
    VectorX<double> Lambda;  // the (requested number of) eigenvalues of the subspace matrix S
    MatrixX<double> Z;  // the (requested number of) eigenvectors of the subspace matrix S

    MatrixX<double> V;  // the subspace of guess vectors in an iterative diagonalization algorithm
    MatrixX<double> VA;  // VA = A * V (implicitly calculated through the matrix-vector product)
    MatrixX<double> X;  // contains the new guesses for the eigenvectors (as a linear combination of the current subspace V)

    MatrixX<double> R;  // the residual vectors
    MatrixX<double> Delta;  // the correction vectors (solutions to the residual equations)


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param A                the matrix whose eigenvalue problem should be solved
     */
    EigenproblemEnvironment(const SquareMatrix<double>& A) :
        A (A)
    {}

    /**
     *  @param matrix_vector_product            a vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix)
     *  @param diagonal                         the diagonal of the matrix whose eigenvalue problem should be solved
     */
    EigenproblemEnvironment(const VectorFunction<double>& matrix_vector_product_function, const VectorX<double>& diagonal, const MatrixX<double>& V) :
        dimension (diagonal.size()),
        matrix_vector_product_function (matrix_vector_product_function),
        diagonal (diagonal),
        V (V),
        VA (MatrixX<double>::Zero(V.rows(), 0))  // the initial environment should have no columns in VA
    {}


    /*
     *  STATIC PUBLIC METHODS
     */

    /**
     *  @param A                the matrix whose eigenvalue problem should be solved
     * 
     *  @return an environment that can be used to solve the dense eigenvalue problem for the given square matrix
     */
    static EigenproblemEnvironment Dense(const SquareMatrix<double>& A) {
        return EigenproblemEnvironment(A);
    }

    /**
     *  @param matrix_vector_product            a vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix)
     *  @param diagonal                         the diagonal of the matrix whose eigenvalue problem should be solved
     *  @param V                                a matrix of initial guess vectors (each column of the matrix is an initial guess vector)
     * 
     *  @return an environment that can be used to solve the eigenvalue problem for the matrix that is represented by the given matrix-vector product
     */
    static EigenproblemEnvironment Iterative(const VectorFunction<double>& matrix_vector_product_function, const VectorX<double>& diagonal, const MatrixX<double>& V) {
        return EigenproblemEnvironment(matrix_vector_product_function, diagonal, V);
    }


    /**
     *  @param A                                the matrix whose eigenvalue problem should be solved
     *  @param V                                a matrix of initial guess vectors (each column of the matrix is an initial guess vector)
     * 
     *  @return an environment that can be used to solve the eigenvalue problem for the matrix that is represented by the given matrix-vector product
     */
    static EigenproblemEnvironment Iterative(const SquareMatrix<double>& A, const MatrixX<double>& V) {

        const auto matrix_vector_product_function = [A](const VectorX<double>& x) { return A * x; };
        return EigenproblemEnvironment::Iterative(matrix_vector_product_function, A.diagonal(), V);  
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param number_of_requested_eigenpairs               the number of eigenpairs you would like to retrieve
     * 
     *  @return the eigenvalues and eigenvectors as a vector of eigenpairs
     */
    std::vector<Eigenpair> eigenpairs(const size_t number_of_requested_eigenpairs = 1) const {

        if (number_of_requested_eigenpairs > eigenvectors.cols()) {
            throw std::invalid_argument("EigenproblemEnvironment::eigenpairs(const size_t): You cannot retrieve that many eigenpairs.");
        }

        std::vector<Eigenpair> eigenpairs {};
        eigenpairs.reserve(number_of_requested_eigenpairs);

        for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
            const auto& eigenvalue = this->eigenvalues(i);
            const auto& eigenvector = this->eigenvectors.col(i);

            eigenpairs.emplace_back(eigenvalue, eigenvector);
        }

        return eigenpairs;
    }
};



}  // namespace GQCP
