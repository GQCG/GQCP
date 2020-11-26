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


#include "Basis/Transformations/RTransformation.hpp"
#include "DensityMatrix/DensityMatrixTraits.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "DensityMatrix/Simple2DM.hpp"


namespace GQCP {


/*
 *  MARK: `Orbital2DM` implementation
 */

/**
 *  A type used to represent a two-electron orbital density matrix, i.e. the summed alpha-alpha, alpha-beta, beta-alpha and beta-beta density matrices.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 */
template <typename _Scalar>
class Orbital2DM:
    public Simple2DM<_Scalar, Orbital2DM<_Scalar>> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `Simple2DM`'s constructors.
    using Simple2DM<Scalar, Orbital2DM<Scalar>>::Simple2DM;
};


/*
 *  MARK: `DensityMatrixTraits`
 */

/**
 *  A type that provides compile-time information on `Orbital2DM` that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
struct DensityMatrixTraits<Orbital2DM<Scalar>> {

    // The type of transformation that is naturally related to an `Orbital2DM`. The only transformations that should be naturally possible for an orbital 2-DM are restricted transformations, thereby assuming that the density matrices for alpha-alpha, alpha-beta, beta-alpha and beta-beta are equal and thus transform similarly.
    using Transformation = RTransformation<Scalar>;

    // The type of the one-electron density matrix that is naturally related to an `Orbital2DM`.
    using OneDM = Orbital1DM<Scalar>;
};


/*
 *  MARK: `BasisTransformableTraits`
 */

/**
 *  A type that provides compile-time information on the abstract interface `BasisTransformable` that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
struct BasisTransformableTraits<Orbital2DM<Scalar>> {

    // The type of transformation that is naturally related to an `Orbital2DM`. The only transformations that should be naturally possible for an orbital 2-DM are restricted transformations, thereby assuming that the density matrices for alpha-alpha, alpha-beta, beta-alpha and beta-beta are equal and thus transform similarly.
    using Transformation = RTransformation<Scalar>;
};


/*
 *  MARK: `JacobiRotatableTraits`
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<Orbital2DM<Scalar>> {

    // The type of Jacobi rotation that is naturally related to a `G2DM`.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
