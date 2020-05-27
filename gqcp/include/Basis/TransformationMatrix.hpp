// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Basis/SpinorBasis/JacobiRotationParameters.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  A class that represents a transformation matrix between two orbital bases. The matrix representation of this transformation matrix is such that a new orbital basis b' is found as
 *      b' = b T ,
 *   in which the current orbitals are collected as elements of a row vector b
 * 
 *  @tparam Scalar            the scalar representation of one of the elements of the transformation matrix
 */
template <typename _Scalar>
class TransformationMatrix: public SquareMatrix<_Scalar> {
public:
    using Scalar = _Scalar;


public:
    /*
     *  CONSTRUCTORS
     */

    using SquareMatrix<Scalar>::SquareMatrix;  // inherit base constructors


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  @param jacobi_rotation_parameters       the parameters that define the Jacobi rotation matrix
     *  @param M                                the dimension of the resulting matrix
     *
     *  @return the corresponding Jacobi rotation matrix. Note that we work with the (cos, sin, -sin, cos) definition
     */
    static TransformationMatrix<Scalar> FromJacobi(const JacobiRotationParameters& jacobi_rotation_parameters, size_t M) {

        double c = std::cos(jacobi_rotation_parameters.angle());
        double s = std::sin(jacobi_rotation_parameters.angle());

        // We'll start the construction with an identity matrix
        TransformationMatrix<Scalar> J = TransformationMatrix<Scalar>::Identity(M, M);

        // And apply the Jacobi rotation as J = I * jacobi_rotation (cfr. B' = B T)
        J.applyOnTheRight(jacobi_rotation_parameters.p(), jacobi_rotation_parameters.q(), Eigen::JacobiRotation<double>(c, s));
        return J;
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  In-place 'transform' this transformation matrix such that the resulting transformation matrix describes this and the other transformation together
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        (*this) = (*this) * T;
    }
};


}  // namespace GQCP
