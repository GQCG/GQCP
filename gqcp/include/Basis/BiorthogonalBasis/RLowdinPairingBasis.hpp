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
#include "Basis/Transformations/RTransformation.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"


namespace GQCP {


/**
 *  MARK: RLowdinPairingBasis implementation
 */


/**
 *  A type used to represent a Lowdin pairing basis, generated from a set of restricted expansion coefficients.
 * 
 *  @tparam _Scalar                 The scalar type used for the expansion coefficients: real or complex.
 */
template <typename _Scalar>
class RLowdinPairingBasis:
    public SimpleLowdinPairingBasis<_Scalar, RLowdinPairingBasis<_Scalar>> {
public:
    // The scalar type used for the biorthogonalized expansion coefficients: real or complex.
    using Scalar = _Scalar;

public:
    /**
     *  MARK: Constructors
     */

    // Inherit `SimpleLowdinPairingBasis`'s constructors.
    using SimpleLowdinPairingBasis<Scalar, RLowdinPairingBasis<Scalar>>::SimpleLowdinPairingBasis;
};


/*
 *  MARK: LowdinPairingBasisTraits
 */

/**
 *  A type that provides compile-time information on Löwdin pairing bases that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar>
struct LowdinPairingBasisTraits<RLowdinPairingBasis<_Scalar>> {
    // The scalar type used to represent a coefficient of the biorthogonalized expansions: real or complex.
    using Scalar = _Scalar;

    // The type of transformation that is naturally related to the `RLowdinPairingBasis`.
    using Transformation = RTransformation<Scalar>;

    // The second-quantized representation of the overlap operator related to the `RLowdinPairingBasis`.
    using SQOverlapOperator = ScalarRSQOneElectronOperator<Scalar>;

    // The type of density matrix naturally associated with the `RLowdinPairingBasis`.
    using DM = Orbital1DM<Scalar>;

    // The type of two-particle density matrix naturally associated with the `RLowdinPairingBasis`.
    using TwoDM = Orbital2DM<Scalar>;
};

}  // namespace GQCP
