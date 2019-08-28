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
#pragma once


#include "Mathematical/SquareMatrix.hpp"
#include "typedefs.hpp"


namespace GQCP {


/**
 *  An extension of the SquareMatrix class with methods of quantum chemical context
 *
 *  @tparam _Scalar      the scalar type
 */
template <typename _Scalar>
class ChemicalMatrix : public SquareMatrix<_Scalar> {
public:
    using Scalar = _Scalar;

    using Base = SquareMatrix<Scalar>;
    using Self = ChemicalMatrix<Scalar>;


public:

    /*
     *  CONSTRUCTORS
     */

    using SquareMatrix<Scalar>::SquareMatrix;  // use base constructors


    /*
     *  GETTERS
     */


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of this matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    size_t dimension() const {
        return this->cols();
    }

    size_t get_K() const { return this->dimension(); };


    /**
     *  In-place transform this chemical matrix to another basis
     * 
     *  @param T                            the transformation matrix
     */
    void basisTransformInPlace(const SquareMatrix<Scalar>& T) {

        (*this) = Self(T.adjoint() * (*this) * T);  // has no aliasing issues (https://eigen.tuxfamily.org/dox/group__TopicAliasing.html)
    }


    /**
     *  In-place transform this chemical matrix to another basis
     * 
     *  @param T                            the transformation matrix
     */
    void basisRotateInPlace(const SquareMatrix<Scalar>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("ChemicalMatrix::basisRotateInPlace(const SquareMatrix<Scalar>&): The given transformation matrix is not unitary.");
        }

        this->basisTransformInPlace(U);
    }


    /**
     *  In-place rotate this chemical matrix using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters
     * 
     *  @note This function should only be available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    void basisRotateInPlace(const JacobiRotationParameters& jacobi_rotation_parameters) {

        const auto p = jacobi_rotation_parameters.get_p();
        const auto q = jacobi_rotation_parameters.get_q();
        const auto angle = jacobi_rotation_parameters.get_angle();

        const double c = std::cos(angle);
        const double s = std::sin(angle);


        // Use Eigen's Jacobi module to apply the Jacobi rotations directly (cfr. T.adjoint() * M * T)
        Eigen::JacobiRotation<double> jacobi (c, s);

        this->applyOnTheLeft(p, q, jacobi.adjoint());
        this->applyOnTheRight(p, q, jacobi);
    }
};


}  // namespace GQCP
