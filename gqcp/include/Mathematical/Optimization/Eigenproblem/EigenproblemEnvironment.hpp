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
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {



/**
 *  An environment used to solve eigenvalue problems for self-adjoint matrices.
 */
class EigenproblemEnvironment {
public:
    VectorFunction<double> matrix_vector_product_function;  // a vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix)

    VectorX<double> diagonal;  // the diagonal of the matrix whose eigenvalue problem should be solved
    std::vector<Eigenpair> eigenpairs;

    SquareMatrix<double> A;  // the matrix whose eigenvalue problem should be solved


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
    EigenproblemEnvironment(const VectorFunction<double>& matrix_vector_product_function, const VectorX<double>& diagonal) :
        matrix_vector_product_function (matrix_vector_product_function),
        diagonal (diagonal)
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
     * 
     *  @return an environment that can be used to solve the eigenvalue problem for the matrix that is represented by the given matrix-vector product
     */
    static EigenproblemEnvironment Iterative(const VectorFunction<double>& matrix_vector_product_function, const VectorX<double>& diagonal) {
        return EigenproblemEnvironment(matrix_vector_product_function, diagonal);
    }
};



}  // namespace GQCP
