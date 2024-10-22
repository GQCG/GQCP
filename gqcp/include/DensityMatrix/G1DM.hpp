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


#include "Basis/Transformations/GTransformation.hpp"
#include "DensityMatrix/DensityMatrixTraits.hpp"
#include "DensityMatrix/Simple1DM.hpp"
#include "DensityMatrix/SpinResolved1DM.hpp"


namespace GQCP {


/*
 *  MARK: G1DM implementation
 */

/**
 *  A type used to represent a one-electron general(ized) density matrix, i.e. the full spinor two-component one-electron density matrix.
 *
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 */
template <typename _Scalar>
class G1DM:
    public Simple1DM<_Scalar, G1DM<_Scalar>> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `Simple1DM`'s constructors.
    using Simple1DM<Scalar, G1DM<Scalar>>::Simple1DM;


    /**
     *  Create a `G1DM` from a `SpinResolved1DM`.
     *
     *  @param D            The spin-resolved 1-DM.
     *
     *  @return A `G1DM` created from a `SpinResolved1DM`.
     */
    static G1DM<Scalar> FromSpinResolved(const SpinResolved1DM<Scalar>& D) {

        // The goal in this named constructor is to build up the general density matrix (2K x 2K) from the spin-resolved (K_sigma x K_sigma) ones.
        const auto K_alpha = D.alpha().numberOfOrbitals();
        const auto K_beta = D.beta().numberOfOrbitals();
        const auto M = K_alpha + K_beta;

        SquareMatrix<Scalar> D_general = SquareMatrix<Scalar>::Zero(M);

        //      alpha |  0
        //        0   | beta
        D_general.topLeftCorner(K_alpha, K_alpha) = D.alpha().matrix();
        D_general.bottomRightCorner(K_beta, K_beta) = D.beta().matrix();

        return G1DM<Scalar> {D_general};
    }
};


/*
 *  MARK: `BasisTransformableTraits`
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<G1DM<Scalar>> {

    // The type of transformation that is naturally related to a `G1DM`.
    using Transformation = GTransformation<Scalar>;
};


/*
 *  MARK: `JacobiRotatableTraits`
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<G1DM<Scalar>> {

    // The type of Jacobi rotation that is naturally related to a `G1DM`.
    using JacobiRotationType = JacobiRotation;
};


/*
 *  MARK: `DensityMatrixTraits`
 */

/**
 *  A type that provides compile-time information on `G1DM` that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
struct DensityMatrixTraits<G1DM<Scalar>> {

    // The type of transformation that is naturally related to a `G1DM`.
    using Transformation = GTransformation<Scalar>;
};


}  // namespace GQCP
