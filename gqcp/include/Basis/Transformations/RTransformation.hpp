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


#include "Basis/Transformations/ROrbitalRotationGenerators.hpp"
#include "Basis/Transformations/SimpleTransformation.hpp"


namespace GQCP {


/*
 *  MARK: RTransformation implementation
 */

/**
 *  A 'restricted' basis transformation, i.e. a spin-orbital basis transformation where the transformation is applied equally to the alpha- and beta-spin-orbitals.
 * 
 *  In general, we adopt the convention outlined in (https://gqcg-res.github.io/knowdes/spinor-transformations.html), where the new orbitals' coefficients can be found in the respective **column** of the related transformation matrix.
 * 
 *  @tparam _Scalar         The scalar type used for a transformation coefficient: real or complex.
 */
template <typename _Scalar>
class RTransformation:
    public SimpleTransformation<_Scalar, RTransformation<_Scalar>> {
public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit SimpleTransformation' constructors.
    using SimpleTransformation<Scalar, RTransformation<Scalar>>::SimpleTransformation;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<RTransformation<Scalar>> {

    // The type of the transformation for which the basis transformation should be defined. A transformation should naturally be transformable with itself.
    using Transformation = RTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<RTransformation<Scalar>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


/*
 *  MARK: OrbitalRotationGeneratorTraits
 */

/**
 *  A type that provides compile-time information on orbital rotation generators that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
struct OrbitalRotationGeneratorTraits<RTransformation<Scalar>> {

    // The type of orbital rotation generators associated with an `RTransformation`.
    using OrbitalRotationGenerators = ROrbitalRotationGenerators<Scalar>;
};


}  // namespace GQCP
