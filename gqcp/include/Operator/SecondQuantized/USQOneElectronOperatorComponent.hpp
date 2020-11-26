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


#include "Basis/MullikenPartitioning/UMullikenPartitioningComponent.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "DensityMatrix/PureSpinResolved2DMComponent.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"
#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Operator/SecondQuantized/SimpleSQOneElectronOperator.hpp"


namespace GQCP {


/**
 *  One of the spin components of an unrestricted one-electron operator, i.e. either the alpha or beta part.
 * 
 *  @tparam _Scalar         The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename _Scalar, typename _Vectorizer>
class USQOneElectronOperatorComponent:
    public SimpleSQOneElectronOperator<_Scalar, _Vectorizer, USQOneElectronOperatorComponent<_Scalar, _Vectorizer>> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SimpleSQOneElectronOperator`'s constructors.
    using SimpleSQOneElectronOperator<Scalar, Vectorizer, USQOneElectronOperatorComponent<Scalar, Vectorizer>>::SimpleSQOneElectronOperator;
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like USQOneElectronOperatorComponent, i.e. with scalar-like access.
template <typename Scalar>
using ScalarUSQOneElectronOperatorComponent = USQOneElectronOperatorComponent<Scalar, ScalarVectorizer>;

// A vector-like USQOneElectronOperatorComponent, i.e. with vector-like access.
template <typename Scalar>
using VectorUSQOneElectronOperatorComponent = USQOneElectronOperatorComponent<Scalar, VectorVectorizer>;

// A matrix-like USQOneElectronOperatorComponent, i.e. with matrix-like access.
template <typename Scalar>
using MatrixUSQOneElectronOperatorComponent = USQOneElectronOperatorComponent<Scalar, MatrixVectorizer>;

// A tensor-like USQOneElectronOperatorComponent, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorUSQOneElectronOperatorComponent = USQOneElectronOperatorComponent<Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information (traits) on `USQOneElectronOperatorComponent` that is otherwise not accessible through a public class alias.
 * 
 *  @tparam Scalar          The scalar type used for a single parameter: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename Scalar, typename Vectorizer>
struct OperatorTraits<USQOneElectronOperatorComponent<Scalar, Vectorizer>> {

    // A type that corresponds to the scalar version of the associated component of an unrestricted one-electron operator type.
    using ScalarOperator = ScalarUSQOneElectronOperatorComponent<Scalar>;

    // The type of transformation that is naturally associated to a component of an unestricted one-electron operator.
    using Transformation = UTransformationComponent<Scalar>;

    // The type of the one-particle density matrix that is naturally associated a component of an unrestricted one-electron operator.
    using OneDM = SpinResolved1DMComponent<Scalar>;

    // The type of the two-particle density matrix that is naturally associated a component of an unrestricted one-electron operator.
    using TwoDM = PureSpinResolved2DMComponent<Scalar>;

    // The type used to encapsulate the Mulliken partitioning scheme.
    using MullikenPartitioning = UMullikenPartitioningComponent<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar, typename Vectorizer>
struct BasisTransformableTraits<USQOneElectronOperatorComponent<Scalar, Vectorizer>> {

    // The type of transformation that is naturally associated to a component of an unrestricted one-electron operator.
    using Transformation = UTransformationComponent<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar, typename Vectorizer>
struct JacobiRotatableTraits<USQOneElectronOperatorComponent<Scalar, Vectorizer>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
