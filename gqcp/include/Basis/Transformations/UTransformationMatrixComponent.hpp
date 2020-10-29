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


#include "Basis/Transformations/SimpleTransformationMatrix.hpp"


namespace GQCP {


/*
 *  MARK: UTransformationMatrixComponent implementation
 */

/**
 *  One of the spin components of an UTransformationMatrix.
 * 
 *  It is specifically designed as one of these spin components, in order to ensuring compile-time correctness. It would be wrong to use either R/GTransformationMatrix as one of the spin components, and it's not possible to use SimpleTransformationMatrix as one of the spin components because it requires a template argument of the type that derives from it.
 * 
 *  @tparam _Scalar         The scalar type used for a transformation coefficient: real or complex.
 */
template <typename _Scalar>
class UTransformationMatrixComponent:
    public SimpleTransformationMatrix<_Scalar, UTransformationMatrixComponent<_Scalar>> {
public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit SimpleTransformationMatrix' constructors.
    using SimpleTransformationMatrix<Scalar, UTransformationMatrixComponent<Scalar>>::SimpleTransformationMatrix;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<UTransformationMatrixComponent<Scalar>> {

    // The type of the transformation matrix for which the basis transformation should be defined. // TODO: Rename "TM" to "TransformationMatrix". A transformation matrix should naturally be transformable with itself.
    using TM = UTransformationMatrixComponent<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<UTransformationMatrixComponent<Scalar>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
