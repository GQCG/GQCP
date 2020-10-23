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


#include "Operator/SecondQuantized/SimpleSQTwoElectronOperator.hpp"


namespace GQCP {


/**
 *  One of the spin components of an unrestricted two-electron operator.
 * 
 *  @tparam _Scalar         The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
 */
template <typename _Scalar, typename _Vectorizer>
class USQTwoElectronOperatorComponent:
    public SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, USQTwoElectronOperator<_Scalar, _Vectorizer>> {
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
    using SimpleSQTwoElectronOperator<Scalar, Vectorizer, USQTwoElectronOperatorComponent<Scalar, Vectorizer>>::SimpleSQTwoElectronOperator;
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like USQTwoElectronOperatorComponent, i.e. with scalar-like access.
template <typename Scalar>
using ScalarUSQTwoElectronOperatorComponent = USQTwoElectronOperatorComponent<Scalar, ScalarVectorizer>;

// A vector-like USQTwoElectronOperatorComponent, i.e. with vector-like access.
template <typename Scalar>
using VectorUSQTwoElectronOperatorComponent = USQTwoElectronOperatorComponent<Scalar, VectorVectorizer>;

// A matrix-like USQTwoElectronOperatorComponent, i.e. with matrix-like access.
template <typename Scalar>
using MatrixUSQTwoElectronOperatorComponent = USQTwoElectronOperatorComponent<Scalar, MatrixVectorizer>;

// A tensor-like USQTwoElectronOperatorComponent, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorUSQTwoElectronOperatorComponent = USQTwoElectronOperatorComponent<Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information (traits) on `USQTwoElectronOperatorComponent` that is otherwise not accessible through a public class alias.
 * 
 *  @tparam Scalar          The scalar type used for a single parameter: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
 */
template <typename Scalar, typename Vectorizer>
struct OperatorTraits<USQTwoElectronOperatorComponent<Scalar, Vectorizer>> {

    // A type that corresponds to the scalar version of the associated component of an unrestricted two-electron operator type.
    using ScalarOperator = ScalarUSQTwoElectronOperatorComponent<Scalar>;

    // The type of transformation matrix that is naturally associated to a component of an unestricted two-electron operator.
    using TM = UTransformationMatrixComponent<Scalar>;

    // The type of the one-particle density matrix that is naturally associated a component of an unestricted two-electron operator.
    using OneDM = SpinResolved1DMComponent<Scalar>;

    // The type of the two-particle density matrix that is naturally associated a component of an unestricted two-electron operator.
    using TwoDM = SpinResolved2DMComponent<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar, typename Vectorizer>
struct BasisTransformableTraits<USQTwoElectronOperatorComponent<Scalar, Vectorizer>> {

    // The type of transformation matrix that is naturally associated to a component of an unestricted two-electron operator.
    using TM = UTransformationMatrixComponent<Scalar>;
};


}  // namespace GQCP
