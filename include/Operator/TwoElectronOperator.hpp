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
#ifndef GQCP_TWOELECTRONOPERATOR_HPP
#define GQCP_TWOELECTRONOPERATOR_HPP


#include <unsupported/Eigen/CXX11/Tensor>

#include "BaseOperator.hpp"


namespace GQCP {



/**
 *  A class that holds the matrix representation of a two-electron operator in an orbital basis
 */
class TwoElectronOperator : public BaseOperator {
private:
    Eigen::Tensor<double, 4> tensor;  // the matrix representation of the two-electron operator


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param tensor
     */
    explicit TwoElectronOperator(const Eigen::Tensor<double, 4>& tensor);


    // OPERATORS
    double operator()(size_t p, size_t q, size_t r, size_t s) const { return this->tensor(p, q, r, s); }


    // GETTERS
    const Eigen::Tensor<double, 4>& get_matrix_representation() const { return this->tensor; }


    // PUBLIC METHODS
    /**
     *  Transform the matrix representation of a two-electron operator using the transformation matrix @param T
     *
     *  Note that the transformation matrix @param T is used as
     *      b' = b T ,
     *  in which the basis functions are collected as elements of a row vector b
     */
    void transform(const Eigen::MatrixXd& T) override;

    /**
     *  Rotate the matrix representation of a two-electron operator using a unitary rotation matrix @param U
     *
     *  Note that the rotation matrix @param U is used as
     *      b' = b U ,
     *  in which the basis functions are collected as elements of a row vector b.
     */
    void rotate(const Eigen::MatrixXd& U) override;

    /**
     *  Rotate the matrix representation of a two-electron operator using the unitary Jacobi rotation matrix U constructed from the @param jacobi_rotation_parameters
     *
     *  Note that
     *      - the rotation matrix @param U is used as
     *          b' = b U ,
     *        in which the basis functions are collected as elements of a row vector b.
     *      - we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters) override;
};



}  // namespace GQCP


#endif  // GQCP_TWOELECTRONOPERATOR_HPP
