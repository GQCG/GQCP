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
#include "miscellaneous.hpp"

#include <iostream>


namespace GQCP {

    
/**
 *  @param jacobi_rotation_parameters       the parameters that define the Jacobi rotation matrix
 *  @param M                                the dimension of the resulting matrix
 *
 *  @return the corresponding Jacobi rotation matrix. Note that we work with the (cos, sin, -sin, cos) definition of the Jacobi rotation matrix
 */
Eigen::MatrixXd jacobiRotationMatrix(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters, size_t M) {

    double c = std::cos(jacobi_rotation_parameters.get_angle());
    double s = std::sin(jacobi_rotation_parameters.get_angle());

    // We'll start the construction with an identity matrix
    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(M, M);

    // And apply the Jacobi rotation as J = I * jacobi_rotation (cfr. B' = B T)
    J.applyOnTheRight(jacobi_rotation_parameters.get_p(), jacobi_rotation_parameters.get_q(), Eigen::JacobiRotation<double> (c, s));
    return J;
}


/**
 *  @param A    the matrix
 *  @param i    row index (starting from 0)
 *  @param j    column index (starting from 0)
 *
 *  @return the i-j minor of the matrix A (i.e. delete the i-th row and j-th column)
 */
Eigen::MatrixXd minor(const Eigen::MatrixXd& A, size_t i, size_t j) {

    // Delete the i-th row
    Eigen::MatrixXd A_i = Eigen::MatrixXd::Zero(A.rows() - 1, A.cols());
    for (size_t i2 = 0; i2 < A.rows(); i2++) {  // loop over A's rows
        if (i2 < i) {
            A_i.row(i2) = A.row(i2);
        } else if (i2 == i) {
            continue;
        } else if (i2 > i) {
            A_i.row(i2-1) = A.row(i2);
        }
    }

    // Delete the j-th column
    Eigen::MatrixXd A_ij = Eigen::MatrixXd::Zero(A.rows() - 1, A.cols() - 1);
    for (size_t j2 = 0; j2 < A.cols(); j2++) {  // loop over A's columns
        if (j2 < j) {
            A_ij.col(j2) = A_i.col(j2);
        } else if (j2 == j) {
            continue;
        } else if (j2 > j) {
            A_ij.col(j2-1) = A_i.col(j2);
        }
    }

    return A_ij;
}

    

/**
 *  @param A        the square matrix
 *
 *  @return the permanent of the given square matrix
 */
double permanent(const Eigen::MatrixXd& A) {

    if (A.rows() != A.cols()) {
        throw std::invalid_argument("The given matrix must be square.");
    }


    // The recursion ends when the given 'matrix' is just a number
    if ((A.rows() == 1) && (A.cols() == 1)) {
        return A(0,0);
    }

    size_t j = 0;  // develop by the first column
    double value = 0.0;
    for (size_t i = 0; i < A.rows(); i++) {
            value += A(i,j) * permanent(minor(A, i,j));
    }

    return value;
}


}  // namespace GQCP
