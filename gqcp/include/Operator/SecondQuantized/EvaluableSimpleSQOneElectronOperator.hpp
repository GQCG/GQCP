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


#include "Mathematical/Functions/Function.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/SQOperatorStorageBase.hpp"


namespace GQCP {


/**
 *  An evaluable, simple second-quantized one-electron operator whose parameters are evaluatable at a point in space.
 *
 *  This type (and its derived types) distinguish(es) itself from `SimpleSQOneElectronOperator`, which is mainly to be used for scalar types that are real or complex, thus enabling conformance to `BasisTransformable` and `JacobiRotatable`. On the contrary, this type does not provide these conformances and transformation formulas should be implemented ad-hoc.
 * 
 *  @tparam _FunctionType                   The type of evaluatable function that is used as a 'matrix element' of this one-electron operator.
 *  @tparam _Vectorizer                     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the density) and vector operators (such as the current density operator).
 *  @tparam _DerivedEvaluableOperator       The type of operator that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _FunctionType, typename _Vectorizer, typename _DerivedEvaluableOperator>
class EvaluableSimpleSQOneElectronOperator:
    public SQOperatorStorageBase<SquareMatrix<_FunctionType>, _Vectorizer, EvaluableSimpleSQOneElectronOperator<_FunctionType, _Vectorizer, _DerivedEvaluableOperator>> {
public:
    // The type of evaluatable function that is used as a matrix element of this one-electron operator.
    using FunctionType = _FunctionType;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the density) and vector operators (such as the current density operator).
    using Vectorizer = _Vectorizer;

    // The type of operator that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedEvaluableOperator = _DerivedEvaluableOperator;

    // The return type of the `operator()`.
    using OutputType = typename FunctionType::OutputType;

    // The input type of the `operator()`.
    using InputType = typename FunctionType::InputType;

    // Allow only `FunctionType` types that actually derive from `Function`.
    static_assert(std::is_base_of<Function<OutputType, InputType>, FunctionType>::value, "ScalarEvaluableRSQOneElectronOperator: FunctionType must inherit from `Function`.");

    // The type that represents the derived operator after its evaluation at a point in space.
    using EvaluatedOperator = typename OperatorTraits<DerivedEvaluableOperator>::EvaluatedOperator;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SQOperatorStorage`'s constructors.
    using SQOperatorStorageBase<SquareMatrix<_FunctionType>, _Vectorizer, EvaluableSimpleSQOneElectronOperator<_FunctionType, _Vectorizer, _DerivedEvaluableOperator>>::SQOperatorStorageBase;


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
    EvaluatedOperator evaluate(const InputType& in) const {

        // Evaluate the operator for every component.
        std::vector<SquareMatrix<OutputType>> F_evaluated_vector {};
        const auto& all_parameters = this->allParameters();
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            SquareMatrix<OutputType> F_evaluated = SquareMatrix<OutputType>::Zero(this->numberOfOrbitals());

            // Evaluate the underlying functions for the given argument.
            for (size_t p = 0; p < this->numberOfOrbitals(); p++) {
                for (size_t q = 0; q < this->numberOfOrbitals(); q++) {
                    F_evaluated(p, q) = all_parameters[i](p, q).operator()(in);  // Evaluate the function at the (p,q)-th element.
                }
            }

            F_evaluated_vector.push_back(F_evaluated);
        }

        const StorageArray<SquareMatrix<OutputType>, Vectorizer> array {F_evaluated_vector, this->vectorizer()};
        return EvaluatedOperator {array};
    }
};


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 * 
 *  @tparam _FunctionType           The type of evaluatable function that is used as a matrix element of this one-electron operator.
 *  @tparam _Vectorizer             The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename _FunctionType, typename _Vectorizer, typename _DerivedEvaluableOperator>
struct OperatorTraits<EvaluableSimpleSQOneElectronOperator<_FunctionType, _Vectorizer, _DerivedEvaluableOperator>> {

    // The type of evaluatable function that is used as a matrix element of the one-electron operator.
    using FunctionType = _FunctionType;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator that derives from `EvaluableSimpleSQOneElectronOperator`, enabling CRTP and compile-time polymorphism.
    using DerivedEvaluableOperator = _DerivedEvaluableOperator;

    // The type of the operator that derives from `EvaluableSimpleSQOneElectronOperator`, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = DerivedEvaluableOperator;
};


}  // namespace GQCP
