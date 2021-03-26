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
 *  @tparam _FunctionType            The type of evaluatable function that is used as a matrix element of this one-electron operator.
 */
template <typename _FunctionType>
class EvaluableScalarRSQOneElectronOperator:
    public SQOperatorStorageBase<SquareMatrix<_FunctionType>, ScalarVectorizer, EvaluableScalarRSQOneElectronOperator<_FunctionType>> {
public:
    // The type of evaluatable function that is used as a matrix element of this one-electron operator.
    using FunctionType = _FunctionType;

    // The return type of the `operator()`.
    using OutputType = typename FunctionType::OutputType;

    // The input type of the `operator()`.
    using InputType = typename FunctionType::InputType;

    // Allow only `FunctionType` types that derive from `Function`.
    static_assert(std::is_base_of<Function<OutputType, InputType>, FunctionType>::value, "EvaluableScalarRSQOneElectronOperator: FunctionType must inherit from `Function`.");

    // The type of 'this'.
    using Self = EvaluableScalarRSQOneElectronOperator<FunctionType>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SQOperatorStorage`'s constructors.
    using SQOperatorStorageBase<SquareMatrix<_FunctionType>, ScalarVectorizer, EvaluableScalarRSQOneElectronOperator<_FunctionType>>::SQOperatorStorageBase;


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
    ScalarRSQOneElectronOperator<OutputType> evaluate(const InputType& x) const {

        // Initialize the results.
        SquareMatrix<OutputType> F_evaluated = SquareMatrix<OutputType>::Zero(this->numberOfOrbitals());

        // Evaluate the underlying scalar functions at the given point.
        for (size_t m = 0; m < this->numberOfOrbitals(); m++) {
            for (size_t n = 0; n < this->numberOfOrbitals(); n++) {
                F_evaluated(m, n) = this->parameters()(m, n).operator()(x);  // Evaluate the function of the (m,n)-th element.
            }
        }

        return ScalarRSQOneElectronOperator<OutputType> {F_evaluated};
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
     *  @note This method is only enabled for EvaluableScalarRSQOneElectronOperator that represent second-quantized electron density operators.
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
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 */
template <typename _FunctionType>
struct OperatorTraits<EvaluableScalarRSQOneElectronOperator<_FunctionType>> {
    // The type of evaluatable function that is used as a matrix element of the one-electron operator.
    using FunctionType = _FunctionType;

    // The scalar type of the evaluatable one-electron operator.
    // using Scalar = typename FunctionType::Scalar;

    // The type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using DerivedOperator = EvaluableScalarRSQOneElectronOperator<FunctionType>;

    // The scalar version of the type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using ScalarOperator = EvaluableScalarRSQOneElectronOperator<FunctionType>;
};


}  // namespace GQCP
