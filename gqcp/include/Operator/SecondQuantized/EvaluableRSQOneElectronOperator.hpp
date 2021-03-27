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
#include "Mathematical/Functions/Function.hpp"
#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQOperatorStorageBase.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A scalar restricted second-quantized one-electron operator whose parameters are evaluatable at a point in space.
 * 
 *  This type distinguishes itself from `RSQOneElectronOperator`, which is mainly to be used for scalar types that are real or complex, thus enabling conformance to `BasisTransformable` and `JacobiRotatable`. On the contrary, this type does not provide these conformances and transformation formulas are implemented ad-hoc.
 * 
 *  @tparam _FunctionType           The type of evaluatable function that is used as a matrix element of this one-electron operator.
 *  @tparam _Vectorizer             The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename _FunctionType, typename _Vectorizer>
class EvaluableRSQOneElectronOperator:
    public SQOperatorStorageBase<SquareMatrix<_FunctionType>, _Vectorizer, EvaluableRSQOneElectronOperator<_FunctionType, _Vectorizer>> {
public:
    // The type of evaluatable function that is used as a matrix element of this one-electron operator.
    using FunctionType = _FunctionType;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The return type of the `operator()`.
    using OutputType = typename FunctionType::OutputType;

    // The input type of the `operator()`.
    using InputType = typename FunctionType::InputType;

    // Allow only `FunctionType` types that derive from `Function`.
    static_assert(std::is_base_of<Function<OutputType, InputType>, FunctionType>::value, "ScalarEvaluableRSQOneElectronOperator: FunctionType must inherit from `Function`.");

    // The type of 'this'.
    using Self = EvaluableRSQOneElectronOperator<FunctionType, Vectorizer>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SQOperatorStorage`'s constructors.
    using SQOperatorStorageBase<SquareMatrix<FunctionType>, Vectorizer, EvaluableRSQOneElectronOperator<FunctionType, Vectorizer>>::SQOperatorStorageBase;


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
    RSQOneElectronOperator<OutputType, Vectorizer> evaluate(const InputType& in) const {

        // Evaluate the operator for every component.
        std::vector<SquareMatrix<OutputType>> F_evaluated_vector {this->numberOfComponents()};
        const auto& all_parameters = this->allParameters();
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            SquareMatrix<OutputType> F_evaluated = SquareMatrix<OutputType>::Zero(this->numberOfOrbitals());

            // Evaluate the underlying functions for the given argument.
            for (size_t p = 0; p < this->numberOfOrbitals(); p++) {
                for (size_t q = 0; q < this->numberOfOrbitals(); q++) {
                    F_evaluated(p, q) = all_parameters[i](p, q).operator()(in);  // Evaluate the function at the (p,q)-th element.
                }
            }

            F_evaluated_vector[i] = F_evaluated;
        }

        return RSQOneElectronOperator<OutputType, Vectorizer> {StorageArray<SquareMatrix<OutputType>, Vectorizer> {F_evaluated_vector, this->vectorizer()}};
    }

    /*
     *  MARK: Calculations
     */

    /**
     *  Evaluate the expectation value of this second-quantized (one-electron) density operator.
     * 
     *  @param D            The 1-DM.
     * 
     *  @return The expectation value of this second-quantized (one-electron) density operator, i.e. the electron density.
     * 
     *  @note This method is only enabled for ScalarEvaluableRSQOneElectronOperator that represent second-quantized electron density operators.
     */
    template <typename S = FunctionType, typename = enable_if_t<std::is_same<S, FunctionProduct<EvaluableLinearCombination<double, EvaluableLinearCombination<double, CartesianGTO>>>>::value>>
    EvaluableLinearCombination<double, FunctionProduct<EvaluableLinearCombination<double, EvaluableLinearCombination<double, CartesianGTO>>>> calculateDensity(const Orbital1DM<double>& D) const {

        using Primitive = CartesianGTO;
        using BasisFunction = EvaluableLinearCombination<double, Primitive>;
        using SpatialOrbital = EvaluableLinearCombination<double, BasisFunction>;
        using SchrodingerDistribution = FunctionProduct<SpatialOrbital>;
        using DensityType = EvaluableLinearCombination<double, SchrodingerDistribution>;

        // Create the density as a linear combination of 'density matrix elements'.
        DensityType density;
        const auto dimension = D.numberOfOrbitals();
        for (size_t p = 0; p < dimension; p++) {
            for (size_t q = 0; q < dimension; q++) {
                const auto coefficient = D.matrix()(p, q);
                const auto function = this->parameters()(p, q);
                density.append(coefficient, function);
            }
        }

        return density;
    }
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

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using DerivedOperator = EvaluableRSQOneElectronOperator<FunctionType, Vectorizer>;

    // The scalar version of the type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using ScalarOperator = ScalarEvaluableRSQOneElectronOperator<FunctionType>;
};


}  // namespace GQCP
