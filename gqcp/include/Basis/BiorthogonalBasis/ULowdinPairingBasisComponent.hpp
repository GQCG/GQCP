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


#include "Basis/BiorthogonalBasis/SimpleLowdinPairingBasis.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"


namespace GQCP {


/**
 *  MARK: ULowdinPairingBasisComponent implementation
 */


/**
 *  A type used to represent a component from a spin resolved Lowdin pairing basis.
 *
 *  @tparam _Scalar                 The scalar type used for the expansion coefficients: real or complex.
 */
template <typename _Scalar>
class ULowdinPairingBasisComponent:
    public SimpleLowdinPairingBasis<_Scalar, ULowdinPairingBasisComponent<_Scalar>> {
public:
    // The scalar type used for the biorthogonalized expansion coefficients: real or complex.
    using Scalar = _Scalar;

public:
    /**
     *  MARK: Constructors
     */

    // Inherit `SimpleLowdinPairingBasis`'s constructors.
    using SimpleLowdinPairingBasis<Scalar, ULowdinPairingBasisComponent<Scalar>>::SimpleLowdinPairingBasis;
};


/*
 *  MARK: LowdinPairingBasisTraits
 */

/**
 *  A type that provides compile-time information on Löwdin pairing bases that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar>
struct LowdinPairingBasisTraits<ULowdinPairingBasisComponent<_Scalar>> {
    // The scalar type used to represent a coefficient of the biorthogonalized expansions: real or complex.
    using Scalar = _Scalar;

    // The type of transformation that is naturally related to the `ULowdinPairingBasisComponent`.
    using Transformation = UTransformationComponent<Scalar>;

    // The second-quantized representation of the overlap operator related to the `ULowdinPairingBasisComponent`.
    using SQOverlapOperator = ScalarUSQOneElectronOperatorComponent<Scalar>;

    // The type of density matrix naturally associated with the `ULowdinPairingBasisComponent`.
    using DM = SpinResolved1DMComponent<Scalar>;

    // The type of density matrix naturally associated with the `ULowdinPairingBasisComponent`.
    using TwoDM = PureSpinResolved2DMComponent<Scalar>;
};

}  // namespace GQCP
