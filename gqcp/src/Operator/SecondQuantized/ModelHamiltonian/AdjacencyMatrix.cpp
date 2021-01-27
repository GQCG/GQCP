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

#include "Operator/SecondQuantized/ModelHamiltonian/AdjacencyMatrix.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  Create an `AdjacencyMatrix` its representation as a `SquareMatrix`.
 * 
 *  @param A        An adjacency matrix, represented as a `SquareMatrix`.
 */
AdjacencyMatrix::AdjacencyMatrix(const SquareMatrix<size_t>& A) :
    A {A} {

    if (!A.isSymmetric()) {
        throw std::invalid_argument("AdjacencyMatrix::AdjacencyMatrix(const SquareMatrix<size_t>&): The given matrix must be symmetric.");
    }

    // Check if the given matrix representation contains either 0s or 1s, and that it has no 1s on the diagonal.
    const std::invalid_argument error {"AdjacencyMatrix::AdjacencyMatrix(const SquareMatrix<size_t>&): The given matrix contains elements that are unexpected for undirected graphs."};
    for (size_t i = 0; i < A.cols(); i++) {
        if (A(i, i) != 0) {
            throw error;
        }

        for (size_t j = 0; j < A.rows(); j++) {
            if (A(i, j) != 0 && A(i, j) != 1) {
                throw error;
            }
        }
    }
}


/*
 *  MARK: Named constructors
 */

/**
 *  @param n            The number of vertices.
 * 
 *  @return An `AdjacencyMatrix` that corresponds to a cyclical undirected graph with n vertices.
 */
AdjacencyMatrix AdjacencyMatrix::Cyclic(const size_t n) {

    // The adjacency matrix of a cyclic graph is the one of a linear graph, where the endpoints are also connected.

    auto A = AdjacencyMatrix::Linear(n);

    A.matrix()(0, n - 1) = 1;
    A.matrix()(n - 1, 0) = 1;

    return A;
}


/**
 *  @param n            The number of vertices.
 * 
 *  @return An `AdjacencyMatrix` that corresponds to a linear undirected graph with n vertices.
 */
AdjacencyMatrix AdjacencyMatrix::Linear(const size_t n) {

    auto A = SquareMatrix<size_t>::Zero(n);

    A.diagonal(1) = VectorX<size_t>::Ones(n - 1);   // The super-diagonal contains 1s.
    A.diagonal(-1) = VectorX<size_t>::Ones(n - 1);  // The sub-diagonal contains 1s.

    return AdjacencyMatrix(A);
}


}  // namespace GQCP
