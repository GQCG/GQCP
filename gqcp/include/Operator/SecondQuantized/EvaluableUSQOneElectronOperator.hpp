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


#include "Operator/SecondQuantized/EvaluableUSQOneElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/**
 *  An unrestricted second-quantized one-electron operator whose parameters are evaluatable at a point in space.
 * 
 *  @tparam _FunctionType                   The type of evaluatable function that is used as a 'matrix element' of this one-electron operator.
 *  @tparam _Vectorizer                     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the density) and vector operators (such as the current density operator).
 */
template <typename _FunctionType, typename _Vectorizer>
class EvaluableUSQOneElectronOperator:
    public SpinResolvedBase<EvaluableUSQOneElectronOperatorComponent<_FunctionType, _Vectorizer>, EvaluableUSQOneElectronOperator<_FunctionType, _Vectorizer>> {
public:
    // The type of evaluatable function that is used as a matrix element of this one-electron operator.
    using FunctionType = _FunctionType;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the density) and vector operators (such as the current density operator).
    using Vectorizer = _Vectorizer;

    // The return type of the `FunctionType::operator()`.
    using OutputType = typename _FunctionType::OutputType;

    // The input type of the `FunctionType::operator()`.
    using InputType = typename _FunctionType::InputType;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<EvaluableUSQOneElectronOperatorComponent<_FunctionType, _Vectorizer>, EvaluableUSQOneElectronOperator<_FunctionType, _Vectorizer>>::SpinResolvedBase;


    /*
     *  MARK: Evaluations
     */

    /**
     *  Evaluate this one-electron operator for the given argument.
     * 
     *  @param in        The argument at which the one-electron operator is to be evaluated.
     *
     *  @return A one-electron operator corresponding to the evaluated functions.
     */
    USQOneElectronOperator<OutputType, Vectorizer> evaluate(const InputType& in) const {

        const auto alpha_evaluation = this->alpha().evaluate(in);
        const auto beta_evaluation = this->beta().evaluate(in);

        return USQOneElectronOperator<OutputType, Vectorizer> {alpha_evaluation, beta_evaluation};
    }
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like EvaluableUSQOneElectronOperator, i.e. with scalar-like access.
template <typename FunctionType>
using ScalarEvaluableUSQOneElectronOperator = EvaluableUSQOneElectronOperator<FunctionType, ScalarVectorizer>;

// A vector-like EvaluableUSQOneElectronOperator, i.e. with vector-like access.
template <typename FunctionType>
using VectorEvaluableUSQOneElectronOperator = EvaluableUSQOneElectronOperator<FunctionType, VectorVectorizer>;

// A matrix-like EvaluableUSQOneElectronOperator, i.e. with matrix-like access.
template <typename FunctionType>
using MatrixEvaluableUSQOneElectronOperator = EvaluableUSQOneElectronOperator<FunctionType, MatrixVectorizer>;

// A tensor-like EvaluableUSQOneElectronOperator, i.e. with tensor-like access.
template <typename FunctionType, size_t N>
using TensorEvaluableUSQOneElectronOperator = EvaluableUSQOneElectronOperator<FunctionType, TensorVectorizer<N>>;


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
struct OperatorTraits<EvaluableUSQOneElectronOperator<_FunctionType, _Vectorizer>> {

    // The type of evaluatable function that is used as a matrix element of the one-electron operator.
    using FunctionType = _FunctionType;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The 1-DM that is naturally associated to the evaluable one-electron operator.
    using OneDM = SpinResolved1DM<OutputType>;
};


}  // namespace GQCP
