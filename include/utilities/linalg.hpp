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
 *  Adds a rank-4 tensor to an other rank-4 tensor starting from given indices
 *
 *  @param T_target         the target rank-4 tensor upon which the addition is performed
 *  @param T                a rank-4 tensor who's value are added to the target
 *
 *  @param i                1st starting index of the starting tensor
 *  @param j                2nd starting index of the starting tensor
 *  @param k                3rd starting index of the starting tensor
 *  @param l                4th starting index of the starting tensor
 */
void tensorBlockAddition(Eigen::Tensor<double, 4>& T_target, const Eigen::Tensor<double, 4>& T, size_t i, size_t j, size_t k, size_t l);


/**
 *  Adds a matrix to a target tensor starting from given indices
 *
 *  @tparam r               indicates which starting index (0,1,2,3) should correspond to the 1st matrix index
 *  @tparam s               indicates which starting index (0,1,2,3) should correspond to the 2nd matrix index
 *
 *  @param T_target         the target rank-4 tensor upon which the addition is performed
 *  @param M                a matrix
 *
 *  @param i                1st starting index of the starting tensor
 *  @param j                2nd starting index of the starting tensor
 *  @param k                3rd starting index of the starting tensor
 *  @param l                4th starting index of the starting tensor
 *
 *
 *  Example:
 *      Given a rank-4 tensor of dimensions (10,10,10,10), and a matrix of dimensions (3,3)
 *      When choosing r and s as 2,1 respectively, and i,j,k,l as 1,2,1,0
 *      The matrix is added to the tensor with 1th index starting at i and 4th index starting at l held fixed at 1 and 0 respectively
 *      While 2nd index starting from j (2) now is mapped to the 2nd matrix (s=1) element and 3rd index starting from k (1) is mapped to 1st (r=2)
 *      meaning that the additions follow : tensor(1,2,1,0) += M(0,0), tensor(1,3,1,0) += M(0,1), ..., tensor(1,3,4,0) += M(3,1), tensor(1,4,4,0) += M(3,2), etc
 *      and for other indices we find : tensor(2,2,1,0) += 0 (no changes, these are obviously not performed)
 */
template<size_t r, size_t s>
void tensorBlockAddition(Eigen::Tensor<double, 4>& T_target, const Eigen::MatrixXd& M, size_t i, size_t j, size_t k, size_t l) {

    // Initialize series of arrays with 1 or 0 values, so that the correct tensor indices given by the template argument correspond to the matrix indices
    size_t ia[4] = {1,0,0,0};
    size_t ja[4] = {0,1,0,0};
    size_t ka[4] = {0,0,1,0};
    size_t la[4] = {0,0,0,1};

    for (size_t x = 0; x < M.rows(); x++) {
        for (size_t y = 0; y < M.cols(); y++) {

            size_t i_effective = i + x * ia[r] + y * ia[s];
            size_t j_effective = j + x * ja[r] + y * ja[s];
            size_t k_effective = k + x * ka[r] + y * ka[s];
            size_t l_effective = l + x * la[r] + y * la[s];

            T_target(i_effective, j_effective, k_effective, l_effective) += M(x,y);
        }
    }
}


/**
 *  Creates a rank-4 tensor from an other rank-4 tensor starting from given indices
 *
 *  @param T                a rank-4 tensor
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
