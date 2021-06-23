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


#include "Basis/Transformations/UTransformationComponent.hpp"
#include "DensityMatrix/Simple2DM.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"


namespace GQCP {


/**
 *  One of the pure (i.e. alpha-alpha or beta-beta) spin components of a spin-resolved 2-DM.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 */
template <typename _Scalar>
class PureSpinResolved2DMComponent:
    public Simple2DM<_Scalar, PureSpinResolved2DMComponent<_Scalar>> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `Simple2DM`'s constructors.
    using Simple2DM<Scalar, PureSpinResolved2DMComponent<Scalar>>::Simple2DM;
};


/*
 *  MARK: `DensityMatrixTraits`
 */

/**
 *  A type that provides compile-time information on `PureSpinResolved2DMComponent` that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
struct DensityMatrixTraits<PureSpinResolved2DMComponent<Scalar>> {

    // The type of transformation that is naturally associated to a component of a spin-resolved 2-DM.
    using Transformation = UTransformationComponent<Scalar>;

    // The type of the one-electron density matrix that is naturally related to an `PureSpinResolved2DMComponent`.
    using OneDM = SpinResolved1DMComponent<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<PureSpinResolved2DMComponent<Scalar>> {

    // The type of transformation that is naturally associated to a component of a spin-resolved 2-DM.
    using Transformation = UTransformationComponent<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<PureSpinResolved2DMComponent<Scalar>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
