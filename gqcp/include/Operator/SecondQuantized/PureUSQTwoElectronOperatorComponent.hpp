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


#include "Basis/Transformations/UTransformationComponent.hpp"
#include "DensityMatrix/PureSpinResolved2DMComponent.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"
#include "Operator/SecondQuantized/SimpleSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"


namespace GQCP {


/**
 *  One of the pure (i.e. alpha-alpha or beta-beta) spin components of an unrestricted two-electron operator.
 * 
 *  @tparam _Scalar         The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
 */
template <typename _Scalar, typename _Vectorizer>
class PureUSQTwoElectronOperatorComponent:
    public SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, PureUSQTwoElectronOperatorComponent<_Scalar, _Vectorizer>> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
    using Vectorizer = _Vectorizer;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SimpleSQTwoElectronOperator`'s constructors.
    using SimpleSQTwoElectronOperator<Scalar, Vectorizer, PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer>>::SimpleSQTwoElectronOperator;
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like PureUSQTwoElectronOperatorComponent, i.e. with scalar-like access.
template <typename Scalar>
using ScalarPureUSQTwoElectronOperatorComponent = PureUSQTwoElectronOperatorComponent<Scalar, ScalarVectorizer>;

// A vector-like PureUSQTwoElectronOperatorComponent, i.e. with vector-like access.
template <typename Scalar>
using VectorPureUSQTwoElectronOperatorComponent = PureUSQTwoElectronOperatorComponent<Scalar, VectorVectorizer>;

// A matrix-like PureUSQTwoElectronOperatorComponent, i.e. with matrix-like access.
template <typename Scalar>
using MatrixPureUSQTwoElectronOperatorComponent = PureUSQTwoElectronOperatorComponent<Scalar, MatrixVectorizer>;

// A tensor-like PureUSQTwoElectronOperatorComponent, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorPureUSQTwoElectronOperatorComponent = PureUSQTwoElectronOperatorComponent<Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information (traits) on `PureUSQTwoElectronOperatorComponent` that is otherwise not accessible through a public class alias.
 * 
 *  @tparam Scalar          The scalar type used for a single parameter: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
 */
template <typename Scalar, typename Vectorizer>
struct OperatorTraits<PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer>> {

    // A type that corresponds to the scalar version of the associated component of an unrestricted two-electron operator type.
    using ScalarOperator = ScalarPureUSQTwoElectronOperatorComponent<Scalar>;

    // The type of one-electron operator that is naturally related to a unrestricted two-electron operator.
    using SQOneElectronOperator = USQOneElectronOperatorComponent<Scalar, Vectorizer>;

    // The type of transformation that is naturally associated to a component of an unrestricted two-electron operator.
    using Transformation = UTransformationComponent<Scalar>;

    // The type of the one-particle density matrix that is naturally associated a component of an unrestricted two-electron operator.
    using OneDM = SpinResolved1DMComponent<Scalar>;

    // The type of the two-particle density matrix that is naturally associated a component of an unrestricted two-electron operator.
    using TwoDM = PureSpinResolved2DMComponent<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar, typename Vectorizer>
struct BasisTransformableTraits<PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer>> {

    // The type of transformation that is naturally associated to a component of an unrestricted two-electron operator.
    using Transformation = UTransformationComponent<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar, typename Vectorizer>
struct JacobiRotatableTraits<PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
