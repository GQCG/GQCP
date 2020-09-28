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


#include "Basis/Transformations/TransformationMatrix.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {


/**
 *  An extension of the SquareMatrix class with methods of quantum chemical context for one-electron integrals.
 *
 *  @tparam _Scalar      the scalar type
 */
template <typename _Scalar>
class QCMatrix: public SquareMatrix<_Scalar> {
public:
    using Scalar = _Scalar;

    using Base = SquareMatrix<Scalar>;
    using Self = QCMatrix<Scalar>;


public:
    /*
     *  CONSTRUCTORS
     */

    using SquareMatrix<Scalar>::SquareMatrix;  // use base constructors


    /*
     *  PUBLIC METHODS
     */

    /**
     *  In-place transform this chemical matrix to another basis
     * 
     *  @param T                            the transformation matrix
     */
    void basisRotate(const TransformationMatrix<Scalar>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("QCMatrix::basisRotate(const TransformationMatrix<Scalar>&): The given transformation matrix is not unitary.");
        }

        this->basisTransform(U);
    }


    /**
     *  In-place rotate this chemical matrix using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters
     * 
     *  @note This function should only be available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    void basisRotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Use Eigen's Jacobi module to apply the Jacobi rotations directly (cfr. T.adjoint() * M * T)
        const auto p = jacobi_rotation_parameters.p();
        const auto q = jacobi_rotation_parameters.q();
        const auto jacobi_rotation = jacobi_rotation_parameters.Eigen();

        this->applyOnTheLeft(p, q, jacobi_rotation.adjoint());
        this->applyOnTheRight(p, q, jacobi_rotation);
    }


    /**
     *  In-place transform this chemical matrix to another basis
     * 
     *  @param T                            the transformation matrix
     */
    void basisTransform(const TransformationMatrix<Scalar>& T) {

        (*this) = Self(T.adjoint() * (*this) * T);  // has no aliasing issues (https://eigen.tuxfamily.org/dox/group__TopicAliasing.html)
    }


    /**
     *  @return the number of orbitals (spinors or spin-orbitals, depending on the context) this quantum chemical matrix is related to
     */
    size_t numberOfOrbitals() const { return this->cols(); }
};


}  // namespace GQCP
