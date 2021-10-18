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
#include "Basis/Transformations/GTransformation.hpp"
#include "DensityMatrix/G1DM.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"


namespace GQCP {


/**
 *  MARK: GLowdinPairingBasis implementation
 */


/**
 *  A type used to represent a Löwdin pairing basis, generated from a set of generalized expansion coefficients.
 * 
 *  @tparam _Scalar                 The scalar type used for the expansion coefficients: real or complex.
 */
template <typename _Scalar>
class GLowdinPairingBasis:
    public SimpleLowdinPairingBasis<_Scalar, GLowdinPairingBasis<_Scalar>> {
public:
    // The scalar type used for the biorthogonalized expansion coefficients: real or complex.
    using Scalar = _Scalar;

public:
    /**
     *  MARK: Constructors
     */

    // Inherit `SimpleLowdinPairingBasis`'s constructors.
    using SimpleLowdinPairingBasis<Scalar, GLowdinPairingBasis<Scalar>>::SimpleLowdinPairingBasis;
};


/*
 *  MARK: SpinorBasisTraits
 */

/**
 *  A type that provides compile-time information on Löwdin pairing bases that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar>
struct LowdinPairingBasisTraits<GLowdinPairingBasis<_Scalar>> {
    // The scalar type used to represent a coefficient of the biorthogonalized expansions: real or complex.
    using Scalar = _Scalar;

    // The type of transformation that is naturally related to a `GLowdinPairingBasis`.
    using Transformation = GTransformation<Scalar>;

    // The second-quantized representation of the overlap operator related to the Löwdin pairing basis.
    using SQOverlapOperator = ScalarGSQOneElectronOperator<Scalar>;

    // The type of matrix naturally associated with a `GLowdinPairingBasis`.
    using Matrix = MatrixX<Scalar>;

    // The type of density matrix naturally associated with a `GLowdinPairingBasis`.
    using DM = G1DM<Scalar>;
};

}  // namespace GQCP
