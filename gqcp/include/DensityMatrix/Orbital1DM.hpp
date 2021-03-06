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
#include "DensityMatrix/Simple1DM.hpp"


namespace GQCP {


/*
 *  MARK: Orbital1DM implementation
 */

/**
 *  A type used to represent a one-electron orbital density matrix, i.e. the summed alpha and beta density matrix.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 * 
 *  @note This type represents the (1/sqrt(2)-scaled) singlet (0,0) 1-DM.
 */
template <typename _Scalar>
class Orbital1DM:
    public Simple1DM<_Scalar, Orbital1DM<_Scalar>> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `Simple1DM`'s constructors.
    using Simple1DM<Scalar, Orbital1DM<Scalar>>::Simple1DM;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<Orbital1DM<Scalar>> {

    // The type of transformation that is naturally related to an `Orbital1DM`.
    using Transformation = RTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<Orbital1DM<Scalar>> {

    // The type of Jacobi rotation that is naturally related to an `Orbital1DM`.
    using JacobiRotationType = JacobiRotation;
};


/*
 *  MARK: DensityMatrixTraits
 */

/**
 *  A type that provides compile-time information on `Orbital1DM` that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
struct DensityMatrixTraits<Orbital1DM<Scalar>> {

    // The type of transformation that is naturally related to an `Orbital1DM`. The only transformations that should be naturally possible for an orbital 1-DM are restricted transformations, thereby assuming that the density matrices for alpha and beta are equal and thus transform similarly.
    using Transformation = RTransformation<Scalar>;
};


}  // namespace GQCP
