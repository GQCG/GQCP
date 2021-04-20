// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Mathematical/Optimization/Eigenproblem/Eigenpair.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  An environment used to solve eigenvalue problems for self-adjoint matrices.
 * 
 *  @tparam _Scalar         The scalar type of the matrix elements: real or complex.
 */
template <typename _Scalar>
class EigenproblemEnvironment {
public:
    // The scalar type of the matrix elements: real or complex.
    using Scalar = _Scalar;


public:
    // A vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix).
    VectorFunction<Scalar> matrix_vector_product_function;


    // The self-adjoint matrix whose eigenvalue problem should be solved.
    SquareMatrix<Scalar> A;

    // The diagonal of the matrix.
    VectorX<Scalar> diagonal;

    // The dimension of the diagonalization problem (the dimension of one eigenvector).
    size_t dimension;


    // The eigenvalues of the matrix A.
    VectorX<double> eigenvalues;

    // The eigenvectors of the matrix A.
    MatrixX<Scalar> eigenvectors;

    // The "subspace matrix": the projection of the matrix A onto the subspace spanned by the vectors in V.
    SquareMatrix<Scalar> S;

    // The (requested number of) eigenvalues of the subspace matrix S.
    VectorX<double> Lambda;

    // The (requested number of) eigenvectors of the subspace matrix S.
    MatrixX<Scalar> Z;


    // The subspace of guess vectors in an iterative diagonalization algorithm.
    MatrixX<Scalar> V;

    // VA = A * V (implicitly calculated through the matrix-vector product).
    MatrixX<Scalar> VA;

    // Contains the new guesses for the eigenvectors (as a linear combination of the current subspace V).
    MatrixX<Scalar> X;


    // The residual vectors.
    MatrixX<Scalar> R;

    // The correction vectors, i.e. the solutions to the residual equations.
    MatrixX<Scalar> Delta;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param A                The matrix whose eigenvalue problem should be solved.
     */
    EigenproblemEnvironment(const SquareMatrix<Scalar>& A) :
        A {A} {}

    /**
     *  @param matrix_vector_product            A vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix).
     *  @param diagonal                         The diagonal of the matrix whose eigenvalue problem should be solved.
     */
    EigenproblemEnvironment(const VectorFunction<Scalar>& matrix_vector_product_function, const VectorX<Scalar>& diagonal, const MatrixX<Scalar>& V) :
        dimension {static_cast<size_t>(diagonal.size())},
        matrix_vector_product_function {matrix_vector_product_function},
        diagonal {diagonal},
        V {V},
        VA {MatrixX<Scalar>::Zero(V.rows(), 0)} {}  // The initial environment should have no columns in VA.


    /*
     *  MARK: Named constructors
     */

    /**
     *  @param A                The matrix whose eigenvalue problem should be solved.
     * 
     *  @return An environment that can be used to solve the dense eigenvalue problem for the given square matrix.
     */
    static EigenproblemEnvironment Dense(const SquareMatrix<Scalar>& A) { return EigenproblemEnvironment(A); }

    /**
     *  @param matrix_vector_product            A vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix).
     *  @param diagonal                         The diagonal of the matrix whose eigenvalue problem should be solved.
     *  @param V                                A matrix of initial guess vectors (each column of the matrix is an initial guess vector).
     * 
     *  @return An environment that can be used to solve the eigenvalue problem for the matrix that is represented by the given matrix-vector product.
     */
    static EigenproblemEnvironment Iterative(const VectorFunction<Scalar>& matrix_vector_product_function, const VectorX<Scalar>& diagonal, const MatrixX<Scalar>& V) { return EigenproblemEnvironment(matrix_vector_product_function, diagonal, V); }

    /**
     *  @param A                                The matrix whose eigenvalue problem should be solved.
     *  @param V                                A matrix of initial guess vectors (each column of the matrix is an initial guess vector).
     * 
     *  @return An environment that can be used to solve the eigenvalue problem for the matrix that is represented by the given matrix-vector product.
     */
    static EigenproblemEnvironment Iterative(const SquareMatrix<Scalar>& A, const MatrixX<Scalar>& V) {

        const auto matrix_vector_product_function = [A](const VectorX<Scalar>& x) { return A * x; };
        return EigenproblemEnvironment::Iterative(matrix_vector_product_function, A.diagonal(), V);
    }


    /*
     *  MARK: Access
     */

    /**
     *  @param number_of_requested_eigenpairs               The number of eigenpairs to retrieve.
     * 
     *  @return The eigenvalues and eigenvectors as a vector of eigenpairs.
     */
    std::vector<Eigenpair<double, Scalar>> eigenpairs(const size_t number_of_requested_eigenpairs = 1) const {

        if (number_of_requested_eigenpairs > this->eigenvectors.cols()) {
            throw std::invalid_argument("EigenproblemEnvironment::eigenpairs(const size_t): You cannot retrieve that many eigenpairs.");
        }

        std::vector<Eigenpair<double, Scalar>> eigenpairs {};
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
