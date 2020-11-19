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


/*
 *  MARK: SpinResolved2DMComponent implementation
 */

/**
 *  One of the spin components of a SpinResolved2DM.
 * 
 *  It is specifically designed as one of these spin components, in order to ensuring compile-time correctness. It would be wrong to use either `Orbital2DM or G2DM as one of the spin components, and it's not possible to use `Simple2DM` as one of the spin components because it requires a template argument of the type that derives from it.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 */
template <typename _Scalar>
class SpinResolved2DMComponent:
    public Simple2DM<_Scalar, SpinResolved2DMComponent<_Scalar>> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `Simple2DM`'s constructors.
    using Simple2DM<Scalar, SpinResolved2DMComponent<Scalar>>::Simple2DM;
};


/*
 *  MARK: DensityMatrixTraits
 */

/**
 *  A type that provides compile-time information on `SpinResolved2DMComponent` that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
struct DensityMatrixTraits<SpinResolved2DMComponent<Scalar>> {
    // The type of transformation that is naturally related to a `SpinResolved2DMComponent`. Since a `SpinResolved2DM` naturally transforms with a `UTransformation`, a `SpinResolved2DMComponent` naturally transforms with a `UTransformationComponent`.
    using Transformation = UTransformationComponent<Scalar>;

    // The type of the one-electron density matrix that is naturally related to a `SpinResolved2DMComponent`. It is the return type of `SpinResolved2DMComponent`'s method `reduce`.
    using OneDM_Placeholder = SpinResolved1DMComponent<Scalar>;
};


}  // namespace GQCP
