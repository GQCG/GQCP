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
#ifndef GQCP_BASEOPERATOR_HPP
#define GQCP_BASEOPERATOR_HPP


#include <Eigen/Dense>

#include "JacobiRotationParameters.hpp"



namespace GQCP {


/**
 *  A base class for the representation of operators in an orbital basis
 */
class BaseOperator {
protected:
    size_t dim;  // dimension of the matrix representation of the operator


public:
    // CONSTRUCTORS
    /**
     * @param dimension     the dimension of the operator, i.e. the number of orbitals
     */
    explicit BaseOperator(size_t dimension);


    // GETTERS
    size_t get_dim() const { return this->dim; }


    // PUBLIC METHODS
    /**
     *  In-place transform the matrix representation of the operator
     *
     *  @param T    the transformation matrix between the old and the new orbital basis, it is used as
     *      b' = b T ,
     *   in which the basis functions are collected as elements of a row vector b
     */
    virtual void transform(const Eigen::MatrixXd& T) = 0;

    /**
     *  In-place rotate the matrix representation of the operator
     *
     *  @param U     the unitary transformation (i.e. rotation) matrix, see transform() for how the transformation matrix between the two bases should be represented
     */
    virtual void rotate(const Eigen::MatrixXd& U) = 0;

    /**
     *  In-place rotate the matrix representation of the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    virtual void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) = 0;
};



}  // namespace GQCP


#endif  // GQCP_BASEOPERATOR_HPP
