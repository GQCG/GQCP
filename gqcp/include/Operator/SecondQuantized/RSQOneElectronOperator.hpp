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


#include "Basis/Transformations/RTransformationMatrix.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Operator/SecondQuantized/SimpleSQOneElectronOperator.hpp"


namespace GQCP {


/**
 *  A restricted one-electron operator, which is suited for expressing non-relativistic (spin-free) one-electron operators.
 * 
 *  @tparam _Scalar         The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename _Scalar, typename _Vectorizer>
class RSQOneElectronOperator:
    public SimpleSQOneElectronOperator<_Scalar, _Vectorizer, RSQOneElectronOperator<_Scalar, _Vectorizer>> {
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
    using SimpleSQOneElectronOperator<_Scalar, _Vectorizer, RSQOneElectronOperator<_Scalar, _Vectorizer>>::SimpleSQOneElectronOperator;
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like RSQOneElectronOperator, i.e. with scalar-like access.
template <typename Scalar>
using ScalarRSQOneElectronOperator = RSQOneElectronOperator<Scalar, ScalarVectorizer>;

// A vector-like RSQOneElectronOperator, i.e. with vector-like access.
template <typename Scalar>
using VectorRSQOneElectronOperator = RSQOneElectronOperator<Scalar, VectorVectorizer>;

// A matrix-like RSQOneElectronOperator, i.e. with matrix-like access.
template <typename Scalar>
using MatrixRSQOneElectronOperator = RSQOneElectronOperator<Scalar, MatrixVectorizer>;

// A tensor-like RSQOneElectronOperator, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorRSQOneElectronOperator = RSQOneElectronOperator<Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information (traits) on `RSQOneElectronOperator` that is otherwise not accessible through a public class alias.
 * 
 *  @tparam _Scalar         The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename _Scalar, typename _Vectorizer>
class OperatorTraits<RSQOneElectronOperator<_Scalar, _Vectorizer>> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The operator whose traits are provided.
    using Operator = RSQOneElectronOperator<Scalar, Vectorizer>;

    // A type that corresponds to the scalar version of the associated restricted one-electron operator type.
    using ScalarOperator = ScalarRSQOneElectronOperator<Scalar>;

    // The type of transformation matrix that is naturally associated to a restricted one-electron operator.
    using TM = RTransformationMatrix<Scalar>;

    // The type of density matrix that is naturally associated to a restricted one-electron operator.
    using OneDM = Orbital1DM<Scalar>;
};


}  // namespace GQCP
