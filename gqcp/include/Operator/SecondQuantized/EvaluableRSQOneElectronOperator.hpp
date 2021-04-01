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


#include "DensityMatrix/Orbital1DM.hpp"
#include "Operator/SecondQuantized/EvaluableSimpleSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"


namespace GQCP {


/**
 *  A restricted second-quantized one-electron operator whose parameters are evaluatable at a point in space.
 * 
 *  @tparam _FunctionType                   The type of evaluatable function that is used as a 'matrix element' of this one-electron operator.
 *  @tparam _Vectorizer                     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the density) and vector operators (such as the current density operator).
 */
template <typename _FunctionType, typename _Vectorizer>
class EvaluableRSQOneElectronOperator:
    public EvaluableSimpleSQOneElectronOperator<_FunctionType, _Vectorizer, EvaluableRSQOneElectronOperator<_FunctionType, _Vectorizer>> {
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
    using EvaluableSimpleSQOneElectronOperator<_FunctionType, _Vectorizer, EvaluableRSQOneElectronOperator<_FunctionType, _Vectorizer>>::EvaluableSimpleSQOneElectronOperator;
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like EvaluableRSQOneElectronOperator, i.e. with scalar-like access.
template <typename FunctionType>
using ScalarEvaluableRSQOneElectronOperator = EvaluableRSQOneElectronOperator<FunctionType, ScalarVectorizer>;

// A vector-like EvaluableRSQOneElectronOperator, i.e. with vector-like access.
template <typename FunctionType>
using VectorEvaluableRSQOneElectronOperator = EvaluableRSQOneElectronOperator<FunctionType, VectorVectorizer>;

// A matrix-like EvaluableRSQOneElectronOperator, i.e. with matrix-like access.
template <typename FunctionType>
using MatrixEvaluableRSQOneElectronOperator = EvaluableRSQOneElectronOperator<FunctionType, MatrixVectorizer>;

// A tensor-like EvaluableRSQOneElectronOperator, i.e. with tensor-like access.
template <typename FunctionType, size_t N>
using TensorEvaluableRSQOneElectronOperator = EvaluableRSQOneElectronOperator<FunctionType, TensorVectorizer<N>>;


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
struct OperatorTraits<EvaluableRSQOneElectronOperator<_FunctionType, _Vectorizer>> {

    // The type of evaluatable function that is used as a matrix element of the one-electron operator.
    using FunctionType = _FunctionType;

    // The return type of the `operator()`.
    using OutputType = typename _FunctionType::OutputType;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using DerivedOperator = EvaluableRSQOneElectronOperator<FunctionType, Vectorizer>;

    // The scalar version of the type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using ScalarOperator = ScalarEvaluableRSQOneElectronOperator<FunctionType>;

    // The type that represents the derived operator after its evaluation at a point in space.
    using EvaluatedOperator = RSQOneElectronOperator<OutputType, Vectorizer>;

    // The 1-DM that is naturally associated to the evaluable one-electron operator.
    using OneDM = Orbital1DM<OutputType>;
};


}  // namespace GQCP
