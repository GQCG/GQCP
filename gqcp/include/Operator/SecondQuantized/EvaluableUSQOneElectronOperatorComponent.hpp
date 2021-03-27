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


#include "Operator/SecondQuantized/EvaluableSimpleSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"


namespace GQCP {


/**
 *  One of the components of an unrestricted second-quantized one-electron operator whose parameters are evaluatable at a point in space.
 * 
 *  @tparam _FunctionType                   The type of evaluatable function that is used as a 'matrix element' of this one-electron operator.
 *  @tparam _Vectorizer                     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the density) and vector operators (such as the current density operator).
 */
template <typename _FunctionType, typename _Vectorizer>
class EvaluableUSQOneElectronOperatorComponent:
    public EvaluableSimpleSQOneElectronOperator<_FunctionType, _Vectorizer, EvaluableUSQOneElectronOperatorComponent<_FunctionType, _Vectorizer>> {
public:
    // The type of evaluatable function that is used as a matrix element of this one-electron operator.
    using FunctionType = _FunctionType;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the density) and vector operators (such as the current density operator).
    using Vectorizer = _Vectorizer;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `EvaluableSimpleSQOneElectronOperator`'s constructors.
    using EvaluableSimpleSQOneElectronOperator<_FunctionType, _Vectorizer, EvaluableUSQOneElectronOperatorComponent<_FunctionType, _Vectorizer>>::EvaluableSimpleSQOneElectronOperator;
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like EvaluableUSQOneElectronOperatorComponent, i.e. with scalar-like access.
template <typename FunctionType>
using ScalarEvaluableUSQOneElectronOperatorComponent = EvaluableUSQOneElectronOperatorComponent<FunctionType, ScalarVectorizer>;

// A vector-like EvaluableUSQOneElectronOperatorComponent, i.e. with vector-like access.
template <typename FunctionType>
using VectorEvaluableUSQOneElectronOperatorComponent = EvaluableUSQOneElectronOperatorComponent<FunctionType, VectorVectorizer>;

// A matrix-like EvaluableUSQOneElectronOperatorComponent, i.e. with matrix-like access.
template <typename FunctionType>
using MatrixEvaluableUSQOneElectronOperatorComponent = EvaluableUSQOneElectronOperatorComponent<FunctionType, MatrixVectorizer>;

// A tensor-like EvaluableUSQOneElectronOperatorComponent, i.e. with tensor-like access.
template <typename FunctionType, size_t N>
using TensorEvaluableUSQOneElectronOperatorComponent = EvaluableUSQOneElectronOperatorComponent<FunctionType, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 * 
 *  @tparam _FunctionType           The type of evaluatable function that is used as a matrix element of this one-electron operator.
 *  @tparam _Vectorizer             The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename _FunctionType, typename _Vectorizer>
struct OperatorTraits<EvaluableUSQOneElectronOperatorComponent<_FunctionType, _Vectorizer>> {

    // The type of evaluatable function that is used as a matrix element of the one-electron operator.
    using FunctionType = _FunctionType;

    // The return type of the `operator()`.
    using OutputType = typename _FunctionType::OutputType;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using DerivedOperator = EvaluableUSQOneElectronOperatorComponent<FunctionType, Vectorizer>;

    // The scalar version of the type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using ScalarOperator = ScalarEvaluableUSQOneElectronOperatorComponent<FunctionType>;

    // The type that represents the derived operator after its evaluation at a point in space.
    using EvaluatedOperator = USQOneElectronOperatorComponent<OutputType, Vectorizer>;
};


}  // namespace GQCP
