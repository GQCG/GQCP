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
#include "DensityMatrix/G1DM.hpp"
#include "DensityMatrix/G2DM.hpp"
#include "Domain/MullikenDomain/GMullikenDomain.hpp"
#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Operator/SecondQuantized/SimpleSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"
#include "QuantumChemical/spinor_tags.hpp"


namespace GQCP {


/**
 *  A general(ized) one-electron operator, which is suited for expressing spin-dependent one-electron operators.
 *
 *  @tparam _Scalar         The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename _Scalar, typename _Vectorizer>
class GSQOneElectronOperator:
    public SimpleSQOneElectronOperator<_Scalar, _Vectorizer, GSQOneElectronOperator<_Scalar, _Vectorizer>> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The spinor tag corresponding to a `GSQOneElectronOperator`.
    using SpinorTag = GeneralSpinorTag;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SimpleSQOneElectronOperator`'s constructors.
    using SimpleSQOneElectronOperator<_Scalar, _Vectorizer, GSQOneElectronOperator<_Scalar, _Vectorizer>>::SimpleSQOneElectronOperator;


    /*
     *  MARK: Named constructors
     */

    /**
     *  Construct a `GSQOneElectronOperator` from a `USQOneElectronOperatorComponent`.
     *
     *  @param f_component          The component of an unrestricted one-electron operator that should be converted.
     */
    static GSQOneElectronOperator<Scalar, Vectorizer> FromUnrestrictedComponent(const USQOneElectronOperatorComponent<Scalar, Vectorizer>& f_component) {

        // We can just wrap the one-electron integrals into the correct class.
        const StorageArray<SquareMatrix<Scalar>, Vectorizer> array {f_component.allParameters(), f_component.vectorizer()};
        return GSQOneElectronOperator<Scalar, Vectorizer> {array};
    }
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like GSQOneElectronOperator, i.e. with scalar-like access.
template <typename Scalar>
using ScalarGSQOneElectronOperator = GSQOneElectronOperator<Scalar, ScalarVectorizer>;

// A vector-like GSQOneElectronOperator, i.e. with vector-like access.
template <typename Scalar>
using VectorGSQOneElectronOperator = GSQOneElectronOperator<Scalar, VectorVectorizer>;

// A matrix-like GSQOneElectronOperator, i.e. with matrix-like access.
template <typename Scalar>
using MatrixGSQOneElectronOperator = GSQOneElectronOperator<Scalar, MatrixVectorizer>;

// A tensor-like GSQOneElectronOperator, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorGSQOneElectronOperator = GSQOneElectronOperator<Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information (traits) on `GSQOneElectronOperator` that is otherwise not accessible through a public class alias.
 *
 *  @tparam Scalar          The scalar type used for a single parameter: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename Scalar, typename Vectorizer>
struct OperatorTraits<GSQOneElectronOperator<Scalar, Vectorizer>> {

    // A type that corresponds to the scalar version of the associated general(ized) one-electron operator type.
    using ScalarOperator = ScalarGSQOneElectronOperator<Scalar>;

    // The type of transformation that is naturally associated to a general(ized) one-electron operator.
    using Transformation = GTransformation<Scalar>;

    // The type of the one-particle density matrix that is naturally associated to a general(ized) one-electron operator.
    using OneDM = G1DM<Scalar>;

    // The type of the two-particle density matrix that is naturally associated a restricted one-electron operator.
    using TwoDM = G2DM<Scalar>;

    // The type used to encapsulate the Mulliken domain.
    using MullikenDomain = GMullikenDomain<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 *
 *  @tparam Scalar          The scalar type used for a single parameter: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename Scalar, typename Vectorizer>
struct BasisTransformableTraits<GSQOneElectronOperator<Scalar, Vectorizer>> {

    // The type of transformation that is naturally associated to a `GSQOneElectronOperator`.
    using Transformation = GTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar, typename Vectorizer>
struct JacobiRotatableTraits<GSQOneElectronOperator<Scalar, Vectorizer>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
