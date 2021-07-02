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
#include "DensityMatrix/G1DM.hpp"
#include "DensityMatrix/Simple2DM.hpp"
#include "DensityMatrix/SpinResolved2DM.hpp"


namespace GQCP {


/*
 *  MARK: G2DM implementation
 */

/**
 *  A type used to represent a two-electron general(ized) density matrix, i.e. the full spinor two-component two-electron density matrix.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 */
template <typename _Scalar>
class G2DM:
    public Simple2DM<_Scalar, G2DM<_Scalar>> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `Simple2DM`'s constructors.
    using Simple2DM<Scalar, G2DM<Scalar>>::Simple2DM;


    /**
     *  Create a `G2DM` from a `SpinResolved2DM`.
     * 
     *  @param d            The spin-resolved 2-DM.
     * 
     *  @return A `G2DM` created from a `SpinResolved2DM`.
     */
    static G2DM<Scalar> FromSpinResolved(const SpinResolved2DM<Scalar>& d) {

        // Prepare some variables.
        const auto& d_aaaa = d.alphaAlpha().tensor();
        const auto& d_aabb = d.alphaBeta().tensor();
        const auto& d_bbaa = d.betaAlpha().tensor();
        const auto& d_bbbb = d.betaBeta().tensor();

        // The goal in this named constructor is to build up the general density matrix from the spin-resolved ones.
        const auto K = d_aaaa.dimension();  // Assume that aaaa, aabb, bbaa and bbbb have the same number of orbitals.
        const auto M = 2 * K;
        SquareRankFourTensor<Scalar> d_generalized = SquareRankFourTensor<Scalar>::Zero(M);

        // Primed indices are indices in the larger representation, normal ones are those in the smaller tensors.
        for (size_t mu_ = 0; mu_ < M; mu_++) {  // mu 'prime'
            const size_t mu = mu_ % K;

            for (size_t nu_ = 0; nu_ < M; nu_++) {  // nu 'prime'
                const size_t nu = nu_ % K;

                for (size_t rho_ = 0; rho_ < M; rho_++) {  // rho 'prime'
                    const size_t rho = rho_ % K;

                    for (size_t lambda_ = 0; lambda_ < M; lambda_++) {  // lambda 'prime'
                        const size_t lambda = lambda_ % K;

                        if ((mu_ < K) && (nu_ < K) && (rho_ < K) && (lambda_ < K)) {
                            d_generalized(mu_, nu_, rho_, lambda_) = d_aaaa(mu, nu, rho, lambda);
                        } else if ((mu_ < K) && (nu_ < K) && (rho_ >= K) && (lambda_ >= K)) {
                            d_generalized(mu_, nu_, rho_, lambda_) = d_aabb(mu, nu, rho, lambda);
                        } else if ((mu_ >= K) && (nu_ >= K) && (rho_ < K) && (lambda_ < K)) {
                            d_generalized(mu_, nu_, rho_, lambda_) = d_bbaa(mu, nu, rho, lambda);
                        } else if ((mu_ >= K) && (nu_ >= K) && (rho_ >= K) && (lambda_ >= K)) {
                            d_generalized(mu_, nu_, rho_, lambda_) = d_bbbb(mu, nu, rho, lambda);
                        }
                    }
                }
            }
        }

        return G2DM<Scalar> {d_generalized};
    }
};


/*
 *  MARK: `DensityMatrixTraits`
 */

/**
 *  A type that provides compile-time information on `G2DM` that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
struct DensityMatrixTraits<G2DM<Scalar>> {

    // The type of transformation that is naturally related to a `G2DM`.
    using Transformation = GTransformation<Scalar>;

    // The type of the one-electron density matrix that is naturally related to a `G2DM`.
    using OneDM = G1DM<Scalar>;
};


/*
 *  MARK: `BasisTransformableTraits`
 */

/**
 *  A type that provides compile-time information on the abstract interface `BasisTransformable` that is otherwise not accessible through a public class alias.
 */
template <typename Scalar>
struct BasisTransformableTraits<G2DM<Scalar>> {

    // The type of transformation that is naturally related to a `G2DM`.
    using Transformation = GTransformation<Scalar>;
};


/*
 *  MARK: `JacobiRotatableTraits`
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<G2DM<Scalar>> {

    // The type of Jacobi rotation that is naturally related to a `G2DM`.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
