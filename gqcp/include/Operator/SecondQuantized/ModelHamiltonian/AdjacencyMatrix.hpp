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


#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  An adjacency matrix for an undirected graph.
 */
class AdjacencyMatrix:
    public SquareMatrix<size_t> {

public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create an `AdjacencyMatrix` its representation as a `SquareMatrix`.
     * 
     *  @param A        An adjacency matrix, represented as a `SquareMatrix`.
     */
    AdjacencyMatrix(const SquareMatrix<size_t>& A);

    /**
     *  The default constructor.
     */
    AdjacencyMatrix();

    /**
     *  A constructor required for compatibility with Pybind11. In its 'Eigen' bindings (eigen.h), it makes a call "Type(fits.rows, fits.cols)". This constructor should be called there.
     */
    AdjacencyMatrix(const size_t cols, const size_t rows);


    /*
     *  MARK: Named constructors
     */

    /**
     *  @param n            The number of vertices.
     * 
     *  @return An `AdjacencyMatrix` that corresponds to a cyclical undirected graph with n vertices.
     */
    static AdjacencyMatrix Cyclic(const size_t n);

    /**
     *  @param n            The number of vertices.
     * 
     *  @return An `AdjacencyMatrix` that corresponds to a linear undirected graph with n vertices.
     */
    static AdjacencyMatrix Linear(const size_t n);
};


}  // namespace GQCP
