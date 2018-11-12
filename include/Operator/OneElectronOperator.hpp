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
#ifndef GQCP_ONEELECTRONOPERATOR_HPP
#define GQCP_ONEELECTRONOPERATOR_HPP


#include <Eigen/Dense>

#include "BaseOperator.hpp"


namespace GQCP {



/**
 *  A class that represents a one-electron operator in an orbital basis
 */
class OneElectronOperator : public BaseOperator {
private:
    Eigen::MatrixXd matrix;  // the matrix representation of the one-electron operator


public:
    // CONSTRUCTORS
    /**
     *  A default constructor setting everything to zero
     */
    OneElectronOperator();

    /**
     *  Constructor based on a given @param matrix
     *
     *  @param matrix   the explicit matrix representation of the one-electron operator
     */
    explicit OneElectronOperator(const Eigen::MatrixXd& matrix);


    // GETTERS
    const Eigen::MatrixXd& get_matrix_representation() const { return this->matrix; }


    // OPERATORS
    /**
     *  @return the matrix element at position (p,q)
     */
    double operator()(size_t p, size_t q) const { return this->matrix(p,q); }

    /**
     *  @param other    the other OneElectronOperator
     *
     *  @return the sum of two OneElectronOperators, i.e. a OneElectronOperator whose matrix representation is the sum of the two matrix representations of the given OneElectronOperators
     */
    GQCP::OneElectronOperator operator+(const GQCP::OneElectronOperator& other);

    /**
     *  @return a OneElectronOperator whose matrix representation is negated
     */
    GQCP::OneElectronOperator operator-();

    /**
     *  @param other    the other OneElectronOperator
     *
     *  @return if the matrix representation of this operator is equal to the matrix representation of the other, within the default tolerance specified by isEqualTo()
     */
    bool operator==(const GQCP::OneElectronOperator& other);


    // PUBLIC METHODS
    /**
     *  @param other        the other OneElectronOperator
     *  @param tolerance    the tolerance for equality of the matrix representations
     *
     *  @return if the matrix representation of this operator is equal to the matrix representation of the other, given a tolerance
     */
    bool isEqualTo(const GQCP::OneElectronOperator& other, double tolerance=1.0e-08) const;

    /**
     *  In-place transform the matrix representation of the one-electron operator
     *
     *  @param T    the transformation matrix between the old and the new orbital basis, it is used as
     *      b' = b T ,
     *   in which the basis functions are collected as elements of a row vector b
     */
    void transform(const Eigen::MatrixXd& T) override;

    /**
     *  In-place rotate the matrix representation of the one-electron operator
     *
     *  @param U     the unitary transformation (i.e. rotation) matrix, see transform() for how the transformation matrix between the two bases should be represented
     */
    void rotate(const Eigen::MatrixXd& U) override;

    /**
     *  In-place rotate the matrix representation of the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    void rotate(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters) override;
};



}  // namespace GQCP



#endif  // GQCP_ONEELECTRONOPERATOR_HPP
