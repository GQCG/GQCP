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
#include "Operator/OneElectronOperator.hpp"

#include <stdexcept>


namespace GQCP {


/*
 *  CONSTRUCTOR
 */

/**
 *  Constructor based on a given @param matrix
 */
OneElectronOperator::OneElectronOperator(const Eigen::MatrixXd& matrix) :
    BaseOperator(matrix.cols()),
    matrix (matrix)
{
    // Check if the one-electron integrals are represented as a square matrix
    if (matrix.cols() != matrix.rows()) {
        throw std::invalid_argument("One-electron integrals have to be represented as a square matrix.");
    }
}



/*
 *  OPERATORS
 */

/**
 *  @return the sum of two OneElectronOperators, i.e. a OneElectronOperator whose matrix representation is the sum
 *  of the two matrix representations of the given OneElectronOperators
 */
GQCP::OneElectronOperator OneElectronOperator::operator+(const GQCP::OneElectronOperator& other) {
    
    return OneElectronOperator(this->matrix + other.matrix);
}


/**
 *  @param other    the other OneElectronOperator
 *
 *  @return if the matrix representation of this operator is equal to the matrix representation of the, within the default tolerance specified by isEqualTo()
 */
bool OneElectronOperator::operator==(const GQCP::OneElectronOperator& other) {
    return this->isEqualTo(other);
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param other        the other OneElectronOperator
 *  @param tolerance    the tolerance for equality of the matrix representations
 *
 *  @return if the matrix representation of this operator is equal to the matrix representation of the other, given a tolerance
 */
bool OneElectronOperator::isEqualTo(const GQCP::OneElectronOperator& other, double tolerance) const {
    return this->matrix.isApprox(other.matrix, tolerance);
}


/**
 *  Transform the matrix representation of a one-electron operator using the transformation matrix @param T
 *
 *  Note that the transformation matrix @param T is used as
 *      b' = b T ,
 *  in which the basis functions are collected as elements of a row vector b
 */
void OneElectronOperator::transform(const Eigen::MatrixXd& T) {
    this->matrix = T.adjoint() * this->matrix * T;
}


/**
 *  Rotate the matrix representation of a one-electron operator using a unitary rotation matrix @param U
 *
 *  Note that the rotation matrix @param U is used as
 *      b' = b U ,
 *  in which the basis functions are collected as elements of a row vector b.
 */
void OneElectronOperator::rotate(const Eigen::MatrixXd& U) {

    // Check if the given matrix is actually unitary
    if (!U.isUnitary(1.0e-12)) {
        throw std::invalid_argument("The given matrix is not unitary.");
    }

    this->transform(U);
}


/**
 *  Rotate the matrix representation of a one-electron operator using the unitary Jacobi rotation matrix U constructed from the @param jacobi_rotation_parameters
 *
 *  Note that
 *      - the rotation matrix @param U is used as
 *          b' = b U ,
 *        in which the basis functions are collected as elements of a row vector b.
 *      - we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
 */
void OneElectronOperator::rotate(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters) {

    auto p = jacobi_rotation_parameters.get_p();
    auto q = jacobi_rotation_parameters.get_q();
    auto angle = jacobi_rotation_parameters.get_angle();

    double c = std::cos(angle);
    double s = std::sin(angle);


    // Use Eigen's Jacobi module to apply the Jacobi rotations directly (cfr. T.adjoint() * h * T)
    Eigen::JacobiRotation<double> jacobi (c, s);

    this->matrix.applyOnTheLeft(p, q, jacobi.adjoint());
    this->matrix.applyOnTheRight(p, q, jacobi);
}



}  // namespace GQCP
