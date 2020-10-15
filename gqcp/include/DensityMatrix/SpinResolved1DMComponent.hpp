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


#include "Basis/Transformations/UTransformationMatrixComponent.hpp"
#include "DensityMatrix/Simple1DM.hpp"


namespace GQCP {


/*
 *  MARK: SpinResolved1DMComponent implementation
 */

/**
 *  One of the spin components of a SpinResolved1DM.
 * 
 *  It is specifically designed as one of these spin components, in order to ensuring compile-time correctness. It would be wrong to use either R/GTransformationMatrix as one of the spin components, and it's not possible to use SimpleTransformationMatrix as one of the spin components because it requires a template argument of the type that derives from it.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 */
template <typename _Scalar>
class SpinResolved1DMComponent:
    public Simple1DM<_Scalar, SpinResolved1DMComponent<_Scalar>> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `Simple1DM`'s constructors.
    using Simple1DM<Scalar, SpinResolved1DMComponent<Scalar>>::Simple1DM;
};


/*
 *  MARK: DensityMatrixTraits
 */

/**
 *  A type that provides compile-time information on `SpinResolved1DMComponent` that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
class DensityMatrixTraits<SpinResolved1DMComponent<Scalar>> {
public:
    // The type of transformation matrix that is naturally related to a `SpinResolved1DMComponent`. Since a `SpinResolved1DM` naturally transforms with a `UTransformationMatrix`, a `SpinResolved1DMComponent` naturally transforms with a `UTransformationMatrixComponent`.
    using TM = UTransformationMatrixComponent<Scalar>;
};


}  // namespace GQCP
