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


#include "Basis/BiorthogonalBasis/GLowdinPairingBasis.hpp"
#include "Basis/NonOrthogonalBasis/SimpleNonOrthogonalStateBasis.hpp"
#include "Basis/Transformations/GTransformation.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/GSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  MARK: GNonOrthogonalStateBasis implementation
 */


/**
 *  A type used to represent a non-orthogonal basis, generated from a set of generalized basis states.
 *
 *  @tparam _Scalar                 The scalar type used for the expansion coefficients of the states: real or complex.
 */
template <typename _Scalar>
class GNonOrthogonalStateBasis:
    public SimpleNonOrthogonalStateBasis<_Scalar, GNonOrthogonalStateBasis<_Scalar>> {
public:
    // The scalar type used for the expansion coefficients of the states: real or complex.
    using Scalar = _Scalar;

public:
    /**
     *  MARK: Constructors
     */

    // Inherit `SimpleNonOrthogonalStateBasis`'s constructors.
    using SimpleNonOrthogonalStateBasis<Scalar, GNonOrthogonalStateBasis<Scalar>>::SimpleNonOrthogonalStateBasis;
};


/*
 *  MARK: NonOrthogonalStateBasisTraits
 */

/**
 *  A type that provides compile-time information on non-orthogonal state bases that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar>
struct NonOrthogonalStateBasisTraits<GNonOrthogonalStateBasis<_Scalar>> {
    // The scalar type used to represent a coefficient of the expansions of the states: real or complex.
    using Scalar = _Scalar;

    // The type of transformation that is naturally related to a `GNonOrthogonalStateBasis`.
    using Transformation = GTransformation<Scalar>;

    // The Biorthogonal Löwdin pairing basis associated with the non-orthogonal state basis..
    using BiorthogonalBasis = GLowdinPairingBasis<Scalar>;
};


/*
 *  MARK: OperatorTraits
 */

/**
 *  A type that provides compile-time information on the operators associated with non-orthogonal state bases that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar>
struct NOSBasisOperatorTraits<GNonOrthogonalStateBasis<_Scalar>> {
    // The scalar type used to represent a coefficient of the expansions of the states: real or complex.
    using Scalar = _Scalar;

    // The second quantized representation of the Hamiltonian that is naturally related to a `GNonOrthogonalStateBasis`.
    using Hamiltonian = GSQHamiltonian<Scalar>;

    // The second-quantized representation of the one-electron operator that is naturally related to a `GNonOrthogonalStateBasis`.
    using OneElectronOperator = ScalarGSQOneElectronOperator<Scalar>;

    // The second-quantized representation of the two-electron operator that is naturally related to a `GNonOrthogonalStateBasis`.
    using TwoElectronOperator = ScalarGSQTwoElectronOperator<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _Scalar>
struct JacobiRotatableTraits<GNonOrthogonalStateBasis<_Scalar>> {

    // The type of Jacobi rotation that is naturally related to a `GNonOrthogonalStateBasis`.
    using JacobiRotationType = JacobiRotation;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _Scalar>
struct BasisTransformableTraits<GNonOrthogonalStateBasis<_Scalar>> {

    // The type of transformation that is naturally related to a `GNonOrthogonalStateBAsis`.
    using Transformation = GTransformation<_Scalar>;
};


}  // namespace GQCP
