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


#include "Mathematical/Functions/ScalarFunction.hpp"
#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
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


    // // The scalar type used for a single parameter: real or complex.
    // using Scalar = _Scalar;

    // // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    // using Vectorizer = _Vectorizer;

    // // The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
    // using DerivedOperator = _DerivedOperator;

    // // The matrix representation of the parameters of (one of the components of) the one-electron operator.
    // using MatrixRepresentation = SquareMatrix<Scalar>;

    // // The type of 'this'.
    // using Self = SimpleSQOneElectronOperator<Scalar, Vectorizer, DerivedOperator>;

    // // The type that corresponds to the scalar version of the derived one-electron operator type.
    // using ScalarDerivedOperator = typename OperatorTraits<DerivedOperator>::ScalarOperator;

    // // The type of transformation matrix that is naturally associated to the derived one-electron operator.
    // using TM = typename OperatorTraits<DerivedOperator>::TM;

    // // The type of density matrix that is naturally associated to the derived one-electron operator.
    // using Derived1DM = typename OperatorTraits<DerivedOperator>::OneDM;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SQOperatorStorage`'s constructors.
    using SQOperatorStorageBase<SquareMatrix<_Evaluatable>, ScalarVectorizer, EvaluatableScalarRSQOneElectronOperator<_Evaluatable>>::SQOperatorStorageBase;


    /*
     *  MARK: Calculations
     */

    /**
     *  Calculate the expectation value of this one-electron operator.
     * 
     *  @param D                The 1-DM that represents the wave function.
     *
     *  @return The expectation value of all components of the one-electron operator.
     */
    // StorageArray<Scalar, Vectorizer> calculateExpectationValue(const Derived1DM& D) const {

    //     if (this->numberOfOrbitals() != D.numberOfOrbitals()) {
    //         throw std::invalid_argument("SimpleSQOneElectronOperator::calculateExpectationValue(const OneDM<Scalar>&): The given 1-DM is not compatible with the one-electron operator.");
    //     }

    //     // Calculate the expectation value for every component of the operator.
    //     const auto& parameters = this->allParameters();
    //     std::vector<Scalar> expectation_values(this->numberOfComponents());  // zero-initialize the vector with a number of elements
    //     for (size_t i = 0; i < this->numberOfComponents(); i++) {
    //         expectation_values[i] = (parameters[i] * D).trace();
    //     }

    //     return StorageArray<Scalar, Vectorizer> {expectation_values, this->array.vectorizer()};
    // }


    /*
     *  MARK: Conforming to BasisTransformable
     */

    /**
     *  Apply the basis transformation and return the resulting one-electron integrals.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed one-electron integrals.
     */
    // DerivedOperator transformed(const TM& transformation_matrix) const override {

    //     // Calculate the basis transformation for every component of the operator.
    //     const auto& parameters = this->allParameters();
    //     auto result = this->allParameters();

    //     for (size_t i = 0; i < this->numberOfComponents(); i++) {
    //         result[i] = transformation_matrix.adjoint() * (parameters[i]) * transformation_matrix;
    //     }

    //     return DerivedOperator {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};
    // }


    // // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    // using BasisTransformable<DerivedOperator, TM>::rotate;

    // // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    // using BasisTransformable<DerivedOperator, TM>::rotated;
};


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 */
template <typename _Evaluatable>
class OperatorTraits<EvaluatableScalarRSQOneElectronOperator<_Evaluatable>> {
public:
    // The type of evaluatable function that is used as a matrix element of the one-electron operator.
    using Evaluatable = _Evaluatable;

    //
    // TODO: Change name to VectorArithmeticScalar?
    using Scalar = typename Evaluatable::Scalar;

    // The origin of this public type name is to let `SQOperatorStorage` implement `VectorSpaceArithmetic` on the final derived type. TODO: Change name to 'final'?
    using DerivedOperator = EvaluatableScalarRSQOneElectronOperator<Evaluatable>;

    // TODO: Change name??
    using ScalarOperator = EvaluatableScalarRSQOneElectronOperator<Evaluatable>;
};


}  // namespace GQCP
