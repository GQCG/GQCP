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


#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/SQOperatorStorage.hpp"


namespace GQCP {


/**
 *  A second-quantized one-electron operator whose parameters are described by a single matrix.
 * 
 *  This class is used as a base class for `RSQOneElectronOperator` and `GSQOneElectronOperator`, since they both admit parameter representations using a single matrix, as opposed to `USQOneElectronOperator`, which uses separate alpha- and beta- matrices. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                 The scalar type used for a single parameter/matrix element: real or complex.
 *  @tparam _Vectorizer             The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 *  @tparam _DerivedOperator        The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
class SimpleSQOneElectronOperator:
    public SQOperatorStorage<SquareMatrix<_Scalar>, _Vectorizer, SimpleSQOneElectronOperator<_Scalar, _Vectorizer, _DerivedOperator>>,
    public BasisTransformable<_DerivedOperator, typename OperatorTraits<_DerivedOperator>::TM>,
    public JacobiRotatable<_DerivedOperator> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = _DerivedOperator;

    // The matrix representation of the parameters of (one of the components of) the one-electron operator.
    using MatrixRepresentation = SquareMatrix<Scalar>;

    // The type of 'this'.
    using Self = SimpleSQOneElectronOperator<Scalar, Vectorizer, DerivedOperator>;

    // The type that corresponds to the scalar version of the derived one-electron operator type.
    using ScalarDerivedOperator = typename OperatorTraits<DerivedOperator>::ScalarOperator;

    // The type of transformation matrix that is naturally associated to the derived one-electron operator.
    using TM = typename OperatorTraits<DerivedOperator>::TM;

    // The type of density matrix that is naturally associated to the derived one-electron operator.
    using Derived1DM = typename OperatorTraits<DerivedOperator>::OneDM;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SQOperatorStorage`'s constructors.
    using SQOperatorStorage<SquareMatrix<_Scalar>, _Vectorizer, SimpleSQOneElectronOperator<_Scalar, _Vectorizer, _DerivedOperator>>::SQOperatorStorage;


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
    StorageArray<Scalar, Vectorizer> calculateExpectationValue(const Derived1DM& D) const {

        if (this->numberOfOrbitals() != D.numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQOneElectronOperator::calculateExpectationValue(const OneDM<Scalar>&): The given 1-DM is not compatible with the one-electron operator.");
        }

        // Calculate the expectation value for every component of the operator.
        const auto& parameters = this->allParameters();
        std::vector<Scalar> expectation_values(this->numberOfComponents());  // zero-initialize the vector with a number of elements
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            expectation_values[i] = (parameters[i] * D).trace();
        }

        return StorageArray<Scalar, Vectorizer> {expectation_values, this->array.vectorizer()};
    }


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
    DerivedOperator transformed(const TM& transformation_matrix) const override {

        // Calculate the basis transformation for every component of the operator.
        const auto& parameters = this->allParameters();
        auto result = this->allParameters();

        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            result[i] = transformation_matrix.adjoint() * (parameters[i]) * transformation_matrix;
        }

        return DerivedOperator {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};
    }


    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<DerivedOperator, TM>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedOperator, TM>::rotated;


    /*
     *  MARK: Conforming to JacobiRotatable
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_parameters        The Jacobi rotation parameters.
     * 
     *  @return The jacobi-transformed object.
     */
    DerivedOperator rotated(const JacobiRotationParameters& jacobi_parameters) const override {

        // Use Eigen's Jacobi module to apply the Jacobi rotations directly (cfr. T.adjoint() * M * T).
        const auto p = jacobi_parameters.p();
        const auto q = jacobi_parameters.q();
        const auto jacobi_rotation = jacobi_parameters.Eigen();

        // Calculate the basis transformation for every component of the operator.
        auto result = this->allParameters();
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            result[i].applyOnTheLeft(p, q, jacobi_rotation.adjoint());
            result[i].applyOnTheRight(p, q, jacobi_rotation);
        }

        return DerivedOperator {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<DerivedOperator>::rotate;
};


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
class OperatorTraits<SimpleSQOneElectronOperator<_Scalar, _Vectorizer, _DerivedOperator>> {
public:
    // The scalar type used for a single parameter/matrix element/integral: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator that derives from `SimpleSQOneElectronOperator`, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = _DerivedOperator;
};

}  // namespace GQCP
