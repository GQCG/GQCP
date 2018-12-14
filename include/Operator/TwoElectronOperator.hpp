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
 *  A class that represents a two-electron operator in an orbital basis
 */
class TwoElectronOperator : public BaseOperator {
private:
    Eigen::Tensor<double, 4> tensor;  // the matrix representation of the two-electron operator


public:
    // CONSTRUCTORS
    /**
     *  @param tensor   the explicit matrix representation of the two-electron operator
     */
    explicit TwoElectronOperator(const Eigen::Tensor<double, 4>& tensor);


    // GETTERS
    const Eigen::Tensor<double, 4>& get_matrix_representation() const { return this->tensor; }


    // OPERATORS
    /**
     *  @return the matrix element at position (p,q,r,s)
     */
    double operator()(size_t p, size_t q, size_t r, size_t s) const { return this->tensor(p, q, r, s); }

    /**
     *  @param other    the other TwoElectronOperator
     *
     *  @return if the matrix representation of this operator is equal to the matrix representation of the other, within the default tolerance specified by isEqualTo()
     */
    bool operator==(const TwoElectronOperator& other) const;


    // PUBLIC METHODS
    /**
     *  @param other        the other TwoElectronOperator
     *  @param tolerance    the tolerance for equality of the matrix representations
     *
     *  @return if the matrix representation of this operator is equal to the matrix representation of the other, given a tolerance
     */
    bool isEqualTo(const TwoElectronOperator& other, double tolerance=1.0e-08) const;

    /**
     *  In-place transform the matrix representation of the two-electron operator
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
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) override;
};



}  // namespace GQCP


#endif  // GQCP_TWOELECTRONOPERATOR_HPP
