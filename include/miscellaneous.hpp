// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef GQCP_MISCELLANEOUS_HPP
#define GQCP_MISCELLANEOUS_HPP


#include <stdlib.h>

#include <Eigen/Dense>

#include "JacobiRotationParameters.hpp"


namespace GQCP {


/**
 *  @param jacobi_rotation_parameters       the parameters that define the Jacobi rotation matrix
 *  @param M                                the dimension of the resulting matrix
 *
 *  @return the corresponding Jacobi rotation matrix. Note that we work with the (cos, sin, -sin, cos) definition of the Jacobi rotation matrix
 */
Eigen::MatrixXd jacobiRotationMatrix(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters, size_t M);

/**
 *  @param A    the matrix
 *  @param i    row index (starting from 0)
 *  @param j    column index (starting from 0)
 *
 *  @return the i-j minor of the matrix A (i.e. delete the i-th row and j-th column)
 */
Eigen::MatrixXd matrixMinor(const Eigen::MatrixXd& A, size_t i, size_t j);

/**
 *  @param A        the square matrix
 *
 *  @return the permanent of the given square matrix using a combinatorial algorithm
 */
double permanent_combinatorial(const Eigen::MatrixXd& A);

/**
 *  @param S    the positive integer to be converted to Gray code
 *
 *  @return the Gray code of the given integer number as a bitset
 */
size_t gray_code(size_t S);

/**
 *  @param A        the square matrix
 *
 *  @return the permanent of the given square matrix using the Ryser algorithm.  Note that this algorithm does not work for dimensions larger than 64 (see https://www.codeproject.com/Articles/21282/%2FArticles%2F21282%2FCompute-Permanent-of-a-Matrix-with-Ryser-s-Algorit)
 */
double permanent_ryser(const Eigen::MatrixXd& A);

/**
 *  @param v            the vector index
 *  @param cols         the number of columns in the matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the row-major major (non-contiguous) index given the corresponding vector index
 */
size_t matrixIndexMajor(size_t v, size_t cols, size_t skipped=0);

/**
 *  @param v            the vector index
 *  @param cols         the number of columns in the matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the row-major minor (contiguous) index given the corresponding vector index
 */
size_t matrixIndexMinor(size_t v, size_t cols, size_t skipped=0);

/**
 *  @param i            the row index
 *  @param j            the column index
 *  @param cols         the number of columns in de matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the vector index given the corresponding row-major matrix indices
 */
size_t vectorIndex(size_t i, size_t j, size_t cols, size_t skipped=0);


}  // namespace GQCP



#endif  // GQCP_MISCELLANEOUS_HPP
