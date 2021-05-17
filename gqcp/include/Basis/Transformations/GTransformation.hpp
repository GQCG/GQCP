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


#include "Basis/Transformations/GOrbitalRotationGenerators.hpp"
#include "Basis/Transformations/SimpleTransformation.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "QuantumChemical/SpinResolved.hpp"


namespace GQCP {


/*
 *  MARK: GTransformation implementation
 */

/**
 *  A 'general' basis transformation, i.e. a general, full-spinor basis transformation where the transformation mixes the alpha- and beta components of the two-component spinors.
 * 
 *  In general, we adopt the convention outlined in (https://gqcg-res.github.io/knowdes/spinor-transformations.html), where the new orbitals' coefficients can be found in the respective **column** of the related transformation matrix.
 * 
 *  @tparam _Scalar         The scalar type used for a transformation coefficient: real or complex.
 */
template <typename _Scalar>
class GTransformation:
    public SimpleTransformation<_Scalar, GTransformation<_Scalar>> {
public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;


private:
    // The number of alpha and beta orbitals.
    SpinResolved<size_t> K;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a `GTransformation` from the transformation matrix that it encapsulates.
     * 
     *  @param T                The transformation matrix that collects the expansion coefficients of the new basis (vectors) in the old basis as columns.
     */
    GTransformation(const SquareMatrix<Scalar>& T) :
        SimpleTransformation<_Scalar, GTransformation<_Scalar>>(T),
        K {SpinResolved<size_t>::FromEqual(T.dimension() / 2)} {}


    /**
     *  Construct a `GTransformation` from the transformation matrix that it encapsulates, where the number of basis functions used for the expansion of the alpha components may differ from the number of basis functions for the beta components.
     * 
     *  @param T                The transformation matrix that collects the expansion coefficients of the new basis (vectors) in the old basis as columns.
     *  @param K_alpha          The number of basis functions that are used for the expansion of the alpha components.
     *  @param K_beta           The number of basis functions that are used for the expansion of the beta components.
     */
    GTransformation(const SquareMatrix<Scalar>& T, const size_t K_alpha, const size_t K_beta) :
        SimpleTransformation<_Scalar, GTransformation<_Scalar>>(T),
        K {K_alpha, K_beta} {

        if (T.dimension() != K_alpha + K_beta) {
            throw std::invalid_argument("GTransformation(const SquareMatrix<Scalar>&, const size_t, const size_t): The given transformation matrix is not compatible with the number of given overbitals.");
        }
    }


    /*
     *  MARK: Named constructors
     */

    /**
     *  Convert an `UTransformation` into its generalized counterpart.
     * 
     *  @param u_transformation             The unrestricted transformation.
     * 
     *  @return The `GTransformation` that corresponds to the `UTransformation`.
     */
    static GTransformation<Scalar> FromUnrestricted(const UTransformation<Scalar>& u_transformation) {

        // The goal in this named constructor is to build up the general coefficient matrix (2K x 2K) from the unrestricted (K_sigma x K_sigma) ones.
        const auto K_alpha = u_transformation.alpha().numberOfOrbitals();
        const auto K_beta = u_transformation.beta().numberOfOrbitals();
        const auto M = K_alpha + K_beta;

        SquareMatrix<Scalar> T_general = SquareMatrix<Scalar>::Zero(M);

        //      alpha |  0
        //        0   | beta
        T_general.topLeftCorner(K_alpha, K_alpha) = u_transformation.alpha().matrix();
        T_general.bottomRightCorner(K_beta, K_beta) = u_transformation.beta().matrix();

        return GTransformation<Scalar> {T_general, K_alpha, K_beta};
    }


    /**
     *  Construct an identity transformation related to different number of alpha and beta atomic spinors.
     */
    static GTransformation<Scalar> Identity(const size_t K_alpha, const size_t K_beta) {

        return GTransformation<Scalar> {SquareMatrix<Scalar>::Identity(K_alpha + K_beta), K_alpha, K_beta};
    }


    /*
     *  MARK: Components
     */

    /**
     *  @return The part of the general transformation matrix that describes the expansion coefficients of the alpha components of the spinors.
     */
    MatrixX<Scalar> alpha() const { return this->matrix().topRows(this->K.alpha()); }

    /**
     *  @return The part of the general transformation matrix that describes the expansion coefficients of the beta components of the spinors.
     */
    MatrixX<Scalar> beta() const { return this->matrix().bottomRows(this->K.beta()); }

    /**
     *  @param sigma            Alpha or beta.
     *
     *  @return The part of the general transformation that describes the spinors of the requested component.
     */
    MatrixX<Scalar> component(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->alpha();
            break;
        }

        case Spin::beta: {
            return this->beta();
            break;
        }
        }
    }
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<GTransformation<Scalar>> {

    // The type of the transformation for which the basis transformation should be defined. A transformation matrix should naturally be transformable with itself.
    using Transformation = GTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<GTransformation<Scalar>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


/*
 *  MARK: OrbitalRotationGeneratorTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `OrbitalRotationGenerators`.
 */
template <typename Scalar>
struct OrbitalRotationGeneratorTraits<GTransformation<Scalar>> {

    // The type of orbital rotation generators associated with an `GTransformation`.
    using OrbitalRotationGenerators = GOrbitalRotationGenerators<Scalar>;
};


}  // namespace GQCP
