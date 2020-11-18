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
#include "Mathematical/Functions/ScalarFunction.hpp"
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
 *  @tparam _Evaluatable            The type of evaluatable function that is used as a matrix element of this one-electron operator.
 */
template <typename _Evaluatable>
class EvaluatableScalarRSQOneElectronOperator:
    public SQOperatorStorageBase<SquareMatrix<_Evaluatable>, ScalarVectorizer, EvaluatableScalarRSQOneElectronOperator<_Evaluatable>> {
public:
    // The type of evaluatable function that is used as a matrix element of this one-electron operator.
    using Evaluatable = _Evaluatable;

    // Allow only `Evaluatable` types that derive from `ScalarFunction`.
    static_assert(std::is_base_of<ScalarFunction<typename Evaluatable::Valued, typename Evaluatable::Scalar, Evaluatable::Cols>, Evaluatable>::value, "EvaluatableScalarRSQOneElectronOperator: Evaluatable must inherit from ScalarFunction.");

    // The type of the scalars of the input vector.
    using Scalar = typename Evaluatable::Scalar;

    // The return type of the scalar function.
    using Valued = typename Evaluatable::Valued;

    // The dimension of the input vector: an integer, or Dynamic representing an unknown number of columns at compile time.
    static constexpr auto Cols = Evaluatable::Cols;

    // The type of 'this'.
    using Self = EvaluatableScalarRSQOneElectronOperator<Evaluatable>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SQOperatorStorage`'s constructors.
    using SQOperatorStorageBase<SquareMatrix<_Evaluatable>, ScalarVectorizer, EvaluatableScalarRSQOneElectronOperator<_Evaluatable>>::SQOperatorStorageBase;


    /*
     *  MARK: Evaluations
     */

    /**
     *  Evaluate this one-electron operator at the given point.
     * 
     *  @param x        The vector/point at which the underlying scalar functions should be evaluated.
     *
     *  @return A one-electron operator corresponding to the evaluated scalar functions.
     */
    ScalarRSQOneElectronOperator<Valued> evaluate(const Vector<Scalar, Cols>& x) const {

        // Initialize the results.
        SquareMatrix<Valued> F_evaluated = SquareMatrix<Valued>::Zero(this->numberOfOrbitals());

        // Evaluate the underlying scalar functions at the given point.
        for (size_t m = 0; m < this->numberOfOrbitals(); m++) {
            for (size_t n = 0; n < this->numberOfOrbitals(); n++) {
                F_evaluated(m, n) = this->parameters()(m, n).operator()(x);  // Evaluate the ScalarFunction of the (m,n)-th element.
            }
        }

        return ScalarRSQOneElectronOperator<Valued> {F_evaluated};
    }

    /*
     *  MARK: Calculations
     */

    /**
     *  Evaluate the expectation value of this second-quantized (one-electron) density operator.
     * 
     *  @param D                the 1-DM
     * 
     *  @return the expectation value of this second-quantized (one-electron) density operator, i.e. the electron density
     * 
     *  @note This method is only enabled for EvaluatableScalarRSQOneElectronOperator that represent second-quantized electron density operators.
     */
    template <typename S = Evaluatable, typename = enable_if_t<std::is_same<S, ScalarFunctionProduct<LinearCombination<double, LinearCombination<double, CartesianGTO>>>>::value>>
    LinearCombination<double, ScalarFunctionProduct<LinearCombination<double, LinearCombination<double, CartesianGTO>>>> calculateDensity(const Orbital1DM<double>& D) const {

        using Primitive = CartesianGTO;
        using BasisFunction = LinearCombination<double, Primitive>;
        using SpatialOrbital = LinearCombination<double, BasisFunction>;
        using SchrodingerDistribution = ScalarFunctionProduct<SpatialOrbital>;
        using DensityType = LinearCombination<double, SchrodingerDistribution>;


        // Create the density as a linear combination of 'density matrix elements'.
        const auto dimension = D.numberOfOrbitals();
        DensityType density;
        for (size_t p = 0; p < dimension; p++) {
            for (size_t q = 0; q < dimension; q++) {
                const auto coefficient = D(p, q);
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
template <typename _Evaluatable>
struct OperatorTraits<EvaluatableScalarRSQOneElectronOperator<_Evaluatable>> {
    // The type of evaluatable function that is used as a matrix element of the one-electron operator.
    using Evaluatable = _Evaluatable;

    // The scalar type of the evaluatable one-electron operator.
    using Scalar = typename Evaluatable::Scalar;

    // The type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using DerivedOperator = EvaluatableScalarRSQOneElectronOperator<Evaluatable>;

    // The scalar version of the type of the operator at the end of the inheritance chain of `SQOperatorStorageBase`.
    using ScalarOperator = EvaluatableScalarRSQOneElectronOperator<Evaluatable>;
};


}  // namespace GQCP
