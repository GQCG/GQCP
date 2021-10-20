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


#include "Basis/BiorthogonalBasis/RLowdinPairingBasis.hpp"
#include "Basis/NonOrthogonalBasis/SimpleNonOrthogonalStateBasis.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  MARK: RNonOrthogonalStateBasis implementation
 */


/**
 *  A type used to represent a non-orthogonal basis, generated from a set of generalized basis states.
 *
 *  @tparam _Scalar                 The scalar type used for the expansion coefficients of the states: real or complex.
 */
template <typename _Scalar>
class RNonOrthogonalStateBasis:
    public SimpleNonOrthogonalStateBasis<_Scalar, RNonOrthogonalStateBasis<_Scalar>> {
public:
    // The scalar type used for the expansion coefficients of the states: real or complex.
    using Scalar = _Scalar;

public:
    /**
     *  MARK: Constructors
     */

    // Inherit `SimpleNonOrthogonalStateBasis`'s constructors.
    using SimpleNonOrthogonalStateBasis<Scalar, RNonOrthogonalStateBasis<Scalar>>::SimpleNonOrthogonalStateBasis;
};


/*
 *  MARK: NonOrthogonalStateBasisTraits
 */

/**
 *  A type that provides compile-time information on non-orthogonal state bases that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar>
struct NonOrthogonalStateBasisTraits<RNonOrthogonalStateBasis<_Scalar>> {
    // The scalar type used to represent a coefficient of the expansions of the states: real or complex.
    using Scalar = _Scalar;

    // The type of transformation that is naturally related to a `RNonOrthogonalStateBasis`.
    using Transformation = RTransformation<Scalar>;

    // The Biorthogonal Löwdin pairing basis associated with the non-orthogonal state basis..
    using BiorthogonalBasis = RLowdinPairingBasis<Scalar>;
};


/*
 *  MARK: NOSBasisOperatorTraits
 */

/**
 *  A type that provides compile-time information on the operators associated with non-orthogonal state bases that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar>
struct NOSBasisOperatorTraits<RNonOrthogonalStateBasis<_Scalar>> {
    // The scalar type used to represent a coefficient of the expansions of the states: real or complex.
    using Scalar = _Scalar;

    // The second quantized representation of the Hamiltonian that is naturally related to a `RNonOrthogonalStateBasis`.
    using Hamiltonian = RSQHamiltonian<Scalar>;

    // The second-quantized representation of the one-electron operator that is naturally related to a `RNonOrthogonalStateBasis`.
    using OneElectronOperator = ScalarRSQOneElectronOperator<Scalar>;

    // The second-quantized representation of the two-electron operator that is naturally related to a `RNonOrthogonalStateBasis`.
    using TwoElectronOperator = ScalarRSQTwoElectronOperator<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _Scalar>
struct JacobiRotatableTraits<RNonOrthogonalStateBasis<_Scalar>> {

    // The type of Jacobi rotation that is naturally related to a `RNonOrthogonalStateBasis`.
    using JacobiRotationType = JacobiRotation;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _Scalar>
struct BasisTransformableTraits<RNonOrthogonalStateBasis<_Scalar>> {

    // The type of transformation that is naturally related to a `RNonOrthogonalStateBAsis`.
    using Transformation = RTransformation<_Scalar>;
};


}  // namespace GQCP
