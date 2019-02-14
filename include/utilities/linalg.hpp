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
#ifndef GQCP_LINALG_HPP
#define GQCP_LINALG_HPP



#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>



namespace GQCP {


/**
 *  @param M            the first tensor
 *  @param T            the other tensor
 *  @param tolerance    the tolerance for element-wise comparison
 *
 *  @return if two rank-4 tensors are approximately equal with respect to a tolerance
 *
 *  This function is implemented because Eigen::Tensor does not have an isApprox yet.
 */
bool areEqual(const Eigen::Tensor<double, 4>& M, const Eigen::Tensor<double, 4>& T, double tolerance);

/**
 *  @param eigenvalues1     the first set of eigenvalues
 *  @param eigenvalues2     the second set of eigenvalues
 *  @param tolerance        the tolerance for comparison
 *
 *  @return if two sets of eigenvalues are equal within a given tolerance
 */
bool areEqualEigenvalues(const Eigen::VectorXd& eigenvalues1, const Eigen::VectorXd& eigenvalues2, double tolerance);

/**
 *  @param eigenvector1     the first eigenvector
 *  @param eigenvector2     the second eigenvector
 *  @param tolerance        the tolerance for comparison
 *
 *  @return if two eigenvectors are equal within a given tolerance
 */
bool areEqualEigenvectors(const Eigen::VectorXd& eigenvector1, const Eigen::VectorXd& eigenvector2, double tolerance);

/**
 *  @param eigenvectors1        the first set of eigenvectors
 *  @param eigenvectors2        the second set of eigenvectors
 *  @param tolerance            the tolerance for comparison
 *
 *  @return if two sets of eigenvectors are equal within a given tolerance
 */
bool areEqualSetsOfEigenvectors(const Eigen::MatrixXd& eigenvectors1, const Eigen::MatrixXd& eigenvectors2, double tolerance);


/**
 *  @param M        the matrix
 *
 *  @return the strictly lower triangular matrix (i.e. without the diagonal elements) as a vector in column-major form
 *
 *          5
 *          1   5       -> (1, 2, 3)
 *          2   3   5
 */
Eigen::VectorXd strictLowerTriangle(const Eigen::MatrixXd& M);


/**
 *  @param a        the lower triangle of the matrix in column major form
 *
 *  @return a matrix in which the lower triangle is filled with the given vector, and the other elements are set to zero
 */
Eigen::MatrixXd fillStrictLowerTriangle(const Eigen::VectorXd& a);


/**
 *  @param T        the rank-4 tensor
 *
 *  @return the reduced form of a rank-4 tensor. The elements of the tensor are put into the matrix such that
 *      M(m,n) = T(i,j,k,l)
 *
 *  in which
 *      m is calculated from i and j in a column-major way
 *      n is calculated from k and l in a column-major way
 */
Eigen::MatrixXd toMatrix(const Eigen::Tensor<double, 4>& T);


/**
 *  @param T        a rank-4 tensor with equal dimensions
 *
 *  @return the strict "lower triangle" as a matrix in column major form
 *
 *  The matrix indices (m,n) come from the tensor indices (i,j,k,l) and are such that:
 *      - m is compounded in a column major way from i and j, with the restriction i>j
 *      - n is compounded in a column major way from k and l, with the restriction k>l
 */
Eigen::MatrixXd strictLowerTriangle(const Eigen::Tensor<double, 4>& T);


/**
 *  Copies a rank-4 tensor into an other rank-4 tensor starting from given indices
 *
 *  @param T_target         a rank-4 tensor
 *  @param T                a smaller rank-4 tensor
 *
 *  @param i                1st starting index of the starting tensor
 *  @param j                2nd starting index of the starting tensor
 *  @param k                3rd starting index of the starting tensor
 *  @param l                4th starting index of the starting tensor
 */
void tensorBlockAddition(Eigen::Tensor<double, 4>& T_target, const Eigen::Tensor<double, 4>& T, size_t i, size_t j, size_t k, size_t l);


/**
 *  Copies a matrix into a target tensor starting from given indices
 *
 *  @tparam r               indicates which starting index (0,1,2,3) should correspond to the 1st matrix index
 *  @tparam s               indicates which starting index (0,1,2,3) should correspond to the 2nd matrix index
 *
 *  @param T_target         a rank-4 tensor
 *  @param M                a matrix
 *
 *  @param i                1st starting index of the starting tensor
 *  @param j                2nd starting index of the starting tensor
 *  @param k                3rd starting index of the starting tensor
 *  @param l                4th starting index of the starting tensor
 */
template<size_t r, size_t s>
void tensorBlockAddition(Eigen::Tensor<double, 4>& T_target, const Eigen::MatrixXd& M, size_t i, size_t j, size_t k, size_t l) {

    size_t ia[4] = {1,0,0,0};
    size_t ja[4] = {0,1,0,0};
    size_t ka[4] = {0,0,1,0};
    size_t la[4] = {0,0,0,1};

    for (size_t x = 0; x < M.rows(); x++) {
        for (size_t y = 0; y < M.cols(); y++) {
            T_target(i + x * ia[r] + y * ia[s], j + x * ja[r] + y * ja[s], k + x * ka[r] + y * ka[s], l + x * la[r] + y * la[s]) += M(x,y);
        }
    }
}


/**
 *  Creates a rank-4 tensor from an other rank-4 tensor starting from given indices
 *
 *  @param T                a smaller rank-4 tensor
 *
 *  @param i                1st starting index of the starting tensor
 *  @param j                2nd starting index of the starting tensor
 *  @param k                3rd starting index of the starting tensor
 *  @param l                4th starting index of the starting tensor
 *
 *  @param desize           early cut-off of index iteration
 */
Eigen::Tensor<double, 4> tensorBlockCreation(const Eigen::Tensor<double, 4>& T, size_t i, size_t j, size_t k, size_t l, size_t desize = 0);

}  // namespace GQCP



#endif  // GQCP_LINALG_HPP
