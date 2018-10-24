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
     *  Constructor based on a given @param dimension
     */
    explicit BaseOperator(size_t dimension);


    // GETTERS
    size_t get_dim() const { return this->dim; }


    // PUBLIC METHODS
    /**
     *  Transform the matrix representation of an operator using the transformation matrix @param T
     *
     *  Note that the transformation matrix @param T is used as
     *      b' = b T ,
     *  in which the basis functions are collected as elements of a row vector b
     */
    virtual void transform(const Eigen::MatrixXd& T) = 0;

    /**
     *  Rotate the matrix representation of an operator using a unitary rotation matrix @param U
     *
     *  Note that the rotation matrix @param U is used as
     *      b' = b U ,
     *  in which the basis functions are collected as elements of a row vector b.
     */
    virtual void rotate(const Eigen::MatrixXd& U) = 0;

    /**
     *  Rotate the matrix representation of an operator using the unitary Jacobi rotation matrix U constructed from the @param jacobi_rotation_parameters
     *
     *  Note that
     *      - the rotation matrix @param U is used as
     *          b' = b U ,
     *        in which the basis functions are collected as elements of a row vector b.
     *      - we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    virtual void rotate(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters) = 0;
};



}  // namespace GQCP


#endif  // GQCP_BASEOPERATOR_HPP
