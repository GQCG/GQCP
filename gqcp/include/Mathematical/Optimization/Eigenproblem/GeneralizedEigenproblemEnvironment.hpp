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
 *  An environment used to solve generalized eigenvalue problems for matrices.
 *
 *  @tparam _Scalar         The scalar type of the matrix elements: real or complex.
 */
template <typename _Scalar>
class GeneralizedEigenproblemEnvironment {
public:
    // The scalar type of the matrix elements: real or complex.
    using Scalar = _Scalar;


public:
    // The matrix whose eigenvalue problem should be solved.
    SquareMatrix<Scalar> A;

    // The transformation matrix needed to solve the generalized eigenvalue problem.
    SquareMatrix<Scalar> S;

    // The eigenvalues of the matrix A.
    VectorX<double> eigenvalues;

    // The eigenvectors of the matrix A.
    MatrixX<Scalar> eigenvectors;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param A                The matrix whose eigenvalue problem should be solved.
     */
    GeneralizedEigenproblemEnvironment(const SquareMatrix<Scalar>& A, const SquareMatrix<Scalar>& S) :
        A {A},
        S {S} {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  @param A                The matrix whose eigenvalue problem should be solved.
     *
     *  @return An environment that can be used to solve the dense eigenvalue problem for the given square matrix.
     */
    static GeneralizedEigenproblemEnvironment Dense(const SquareMatrix<Scalar>& A, const SquareMatrix<Scalar>& S) { return GeneralizedEigenproblemEnvironment(A, S); }


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
            throw std::invalid_argument("GeneralEigenproblemEnvironment::eigenpairs(const size_t): You cannot retrieve that many eigenpairs.");
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
