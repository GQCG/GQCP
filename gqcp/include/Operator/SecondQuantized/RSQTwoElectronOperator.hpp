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
#include "DensityMatrix/Orbital1DM.hpp"
#include "DensityMatrix/Orbital2DM.hpp"
#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Operator/SecondQuantized/MixedUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/PureUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SimpleSQTwoElectronOperator.hpp"
#include "QuantumChemical/spinor_tags.hpp"


namespace GQCP {


/**
 *  A restricted two-electron operator, which is suited for expressing non-relativistic (spin-free) two-electron operators.
 * 
 *  @tparam _Scalar                 The scalar type used for a single parameter/matrix element: real or complex.
 *  @tparam _Vectorizer             The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
 */
template <typename _Scalar, typename _Vectorizer>
class RSQTwoElectronOperator:
    public SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, RSQTwoElectronOperator<_Scalar, _Vectorizer>> {
public:
    // The scalar type used for a single parameter/matrix element: real or complex.
    using Scalar = _Scalar;

    //The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
    using Vectorizer = _Vectorizer;

    // The spinor tag corresponding to an `RSQTwoElectronOperator`.
    using SpinorTag = RestrictedSpinOrbitalTag;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SimpleSQOneElectronOperator`'s constructors.
    using SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, RSQTwoElectronOperator<_Scalar, _Vectorizer>>::SimpleSQTwoElectronOperator;


    /*
     *  MARK: Conversions to spin components
     */

    /**
     *  @return The alpha-alpha-component of this restricted two-electron operator.
     */
    PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer> alphaAlpha() const {

        // Since g_pqrs = g_{p alpha, q alpha, r alpha, s alpha}, we can just wrap the two-electron integrals into the correct class.
        const StorageArray<SquareRankFourTensor<Scalar>, Vectorizer> array {this->allParameters(), this->vectorizer()};
        return PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer> {array};
    }


    /**
     *  @return The alpha-beta-component of this restricted two-electron operator.
     */
    MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer> alphaBeta() const {

        // Since g_pqrs = g_{p alpha, q alpha, r beta, s beta}, we can just wrap the two-electron integrals into the correct class.
        const StorageArray<SquareRankFourTensor<Scalar>, Vectorizer> array {this->allParameters(), this->vectorizer()};
        return MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer> {array};
    }


    /**
     *  @return The alpha-beta-component of this restricted two-electron operator.
     */
    MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer> betaAlpha() const {

        // The alpha-beta- and beta-alpha integrals are equal for a restricted operator.
        return this->alphaBeta();
    }


    /*
     *  @return The beta-beta-component of this restricted two-electron operator.
     */
    PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer> betaBeta() const {

        // The alpha-alpha- and beta-beta integrals are equal for a restricted operator.
        return this->alphaAlpha();
    }
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like RSQTwoElectronOperator, i.e. with scalar-like access.
template <typename Scalar>
using ScalarRSQTwoElectronOperator = RSQTwoElectronOperator<Scalar, ScalarVectorizer>;

// A vector-like RSQTwoElectronOperator, i.e. with vector-like access.
template <typename Scalar>
using VectorRSQTwoElectronOperator = RSQTwoElectronOperator<Scalar, VectorVectorizer>;

// A matrix-like RSQTwoElectronOperator, i.e. with matrix-like access.
template <typename Scalar>
using MatrixRSQTwoElectronOperator = RSQTwoElectronOperator<Scalar, MatrixVectorizer>;

// A tensor-like RSQTwoElectronOperator, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorRSQTwoElectronOperator = RSQTwoElectronOperator<Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information (traits) on `RSQTwoElectronOperator` that is otherwise not accessible through a public class alias.
 * 
 *  @tparam Scalar          The scalar type used for a single parameter/matrix element: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
 */
template <typename Scalar, typename Vectorizer>
struct OperatorTraits<RSQTwoElectronOperator<Scalar, Vectorizer>> {

    // A type that corresponds to the scalar version of the associated restricted two-electron operator type.
    using ScalarOperator = ScalarRSQTwoElectronOperator<Scalar>;

    // The type of one-electron operator that is naturally related to a restricted two-electron operator.
    using SQOneElectronOperator = RSQOneElectronOperator<Scalar, Vectorizer>;

    // The type of transformation that is naturally associated to an `RSQTwoElectronOperator`.
    using Transformation = RTransformation<Scalar>;

    // The type of density matrix that is naturally associated to a restricted two-electron operator.
    using OneDM = Orbital1DM<Scalar>;

    // The type of density matrix that is naturally associated to a restricted two-electron operator.
    using TwoDM = Orbital2DM<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 * 
 *  @tparam Scalar          The scalar type used for a single parameter/matrix element: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
 */
template <typename Scalar, typename Vectorizer>
struct BasisTransformableTraits<RSQTwoElectronOperator<Scalar, Vectorizer>> {

    // The type of transformation that is naturally associated to an `RSQTwoElectronOperator`.
    using Transformation = RTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar, typename Vectorizer>
struct JacobiRotatableTraits<RSQTwoElectronOperator<Scalar, Vectorizer>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
