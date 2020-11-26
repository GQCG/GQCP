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
    public BasisTransformable<_DerivedOperator>,
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

    // The type of transformation that is naturally associated to the derived one-electron operator.
    using Transformation = typename OperatorTraits<DerivedOperator>::Transformation;

    // The type of the one-particle density matrix that is naturally associated to the derived one-electron operator.
    using Derived1DM = typename OperatorTraits<DerivedOperator>::OneDM;

    // The type of the two-particle density matrix that is naturally associated to the derived one-electron operator.
    using Derived2DM = typename OperatorTraits<DerivedOperator>::TwoDM;

    // The type used to encapsulate the Mulliken partitioning scheme.
    using MullikenPartitioning = typename OperatorTraits<DerivedOperator>::MullikenPartitioning;


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
     *  @param D                The 1-DM (that represents the wave function).
     *
     *  @return The expectation value of all components of the one-electron operator.
     */
    StorageArray<Scalar, Vectorizer> calculateExpectationValue(const Derived1DM& D) const {

        if (this->numberOfOrbitals() != D.numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQOneElectronOperator::calculateExpectationValue(const Derived1DM<Scalar>&): The given 1-DM's dimension is not compatible with the one-electron operator.");
        }

        // Calculate the expectation value for every component of the operator.
        const auto& parameters = this->allParameters();
        std::vector<Scalar> expectation_values(this->numberOfComponents());  // Zero-initialize the vector with a number of elements.
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            expectation_values[i] = (parameters[i].transpose() * D).trace();
        }

        return StorageArray<Scalar, Vectorizer> {expectation_values, this->array.vectorizer()};
    }


    /**
     *  Calculate the Fockian matrix for (each of the components of) this one-electron operator.
     * 
     *  @param D      The 1-DM (or the response 1-DM for made-variational wave function models).
     *  @param d      The 2-DM (or the response 2-DM for made-variational wave function models).
     *
     *  @return The Fockian matrix.
     * 
     *  @note This method is only enabled in the real case.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, StorageArray<SquareMatrix<double>, Vectorizer>> calculateFockianMatrix(const Derived1DM& D, const Derived2DM& d) const {

        if (D.numberOfOrbitals() != this->numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQOneElectronOperator::calculateFockianMatrix(const Derived1DM&, const Derived2DM&): The 1-DM's dimensions are not compatible with this one-electron operator.");
        }

        if (d.numberOfOrbitals() != this->numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQOneElectronOperator::calculateFockianMatrix(const Derived1DM&, const Derived2DM&): The 2-DM's dimensions are not compatible with this one-electron operator.");
        }

        const auto& parameters = this->allParameters();                           // The parameters of the one-electron operator, as a vector.
        std::vector<SquareMatrix<double>> F_vector {this->numberOfComponents()};  // The resulting vector of Fockian matrices.

        // A KISS implementation of the calculation of the Fockian matrix.
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            const auto& f_i = parameters[i];  // The matrix representation of the parameters of the i-th component.

            // Calculate the Fockian matrix for every component and add it to the vector.
            SquareMatrix<double> F_i = SquareMatrix<double>::Zero(this->numberOfOrbitals());  // The Fockian matrix of the i-th component.
            for (size_t p = 0; p < this->numberOfOrbitals(); p++) {
                for (size_t q = 0; q < this->numberOfOrbitals(); q++) {

                    for (size_t r = 0; r < this->numberOfOrbitals(); r++) {
                        F_i(p, q) += f_i(q, r) * 0.5 * (D(p, r) + D(r, p));  // Include a factor 1/2 to accommodate for response density matrices.
                    }
                }
            }  // F_i elements loop
            F_vector[i] = F_i;
        }

        return StorageArray<SquareMatrix<double>, Vectorizer> {F_vector, this->array.vectorizer()};
    }


    /**
     *  Calculate the super-Fockian matrix for (each of the components of) this one-electron operator.
     * 
     *  @param D      The 1-DM (or the response 1-DM for made-variational wave function models).
     *  @param d      The 2-DM (or the response 2-DM for made-variational wave function models).
     *
     *  @return The super-Fockian matrix.
     * 
     *  @note This method is only enabled in the real case.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, StorageArray<SquareRankFourTensor<double>, Vectorizer>> calculateSuperFockianMatrix(const Derived1DM& D, const Derived2DM& d) const {

        if (D.numberOfOrbitals() != this->numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQOneElectronOperator::calculateFockianMatrix(const Derived1DM&, const Derived2DM&): The given 1-DM's dimensions are not compatible with this one-electron operator.");
        }

        if (d.numberOfOrbitals() != this->numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQOneElectronOperator::calculateFockianMatrix(const Derived1DM&, const Derived2DM&): The given 2-DM's dimensions are not compatible with this one-electron operator.");
        }


        const auto& parameters = this->allParameters();                                   // The parameters of the one-electron operator, as a vector.
        std::vector<SquareRankFourTensor<double>> G_vector {this->numberOfComponents()};  // The resulting vector of super-Fockian matrices.
        const auto F_vector = this->calculateFockianMatrix(D, d).elements();              // The Fockian matrices are necessary in the calculation of the super-Fockian matrices.

        // A KISS implementation of the calculation of the super-Fockian matrix.
        for (size_t i = 0; i < this->numberOfComponents(); i++) {

            const auto& f_i = parameters[i];  // The matrix representation of the parameters of the i-th component.
            const auto& F_i = F_vector[i];    // The Fockian matrix of the i-th component.

            // Calculate the super-Fockian matrix for every component and add it to the array. Add factors 1/2 to accommodate for response density matrices.
            SquareRankFourTensor<double> G_i {this->numberOfOrbitals()};
            G_i.setZero();
            for (size_t p = 0; p < this->numberOfOrbitals(); p++) {
                for (size_t q = 0; q < this->numberOfOrbitals(); q++) {
                    for (size_t r = 0; r < this->numberOfOrbitals(); r++) {
                        for (size_t s = 0; s < this->numberOfOrbitals(); s++) {

                            if (q == r) {
                                G_i(p, q, r, s) += F_i(p, s);
                            }

                            G_i(p, q, r, s) -= f_i(s, p) * 0.5 * (D(r, q) + D(q, r));
                        }
                    }
                }
            }  // G_i elements loop
            G_vector[i] = G_i;
        }

        return StorageArray<SquareRankFourTensor<double>, Vectorizer> {G_vector, this->array.vectorizer()};
    }


    /*
     *  MARK: Conforming to BasisTransformable
     */

    /**
     *  Apply the basis transformation and return the resulting one-electron integrals.
     * 
     *  @param T            The basis transformation.
     * 
     *  @return The basis-transformed one-electron integrals.
     */
    DerivedOperator transformed(const Transformation& T) const override {

        // Calculate the basis transformation for every component of the operator.
        const auto& parameters = this->allParameters();
        auto result = this->allParameters();

        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            const auto& f_i = parameters[i];

            result[i] = T.matrix().adjoint() * f_i * T.matrix();
        }

        return DerivedOperator {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};
    }


    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<DerivedOperator>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedOperator>::rotated;


    /*
     *  MARK: Conforming to JacobiRotatable
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_rotation          The Jacobi rotation.
     * 
     *  @return The jacobi-transformed object.
     */
    DerivedOperator rotated(const JacobiRotation& jacobi_rotation) const override {

        // Use Eigen's Jacobi module to apply the Jacobi rotations directly (cfr. T.adjoint() * M * T).
        const auto p = jacobi_rotation.p();
        const auto q = jacobi_rotation.q();
        const auto jacobi_rotation_eigen = jacobi_rotation.Eigen();

        // Calculate the basis transformation for every component of the operator.
        auto result = this->allParameters();
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            result[i].applyOnTheLeft(p, q, jacobi_rotation_eigen.adjoint());
            result[i].applyOnTheRight(p, q, jacobi_rotation_eigen);
        }

        return DerivedOperator {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<DerivedOperator>::rotate;


    /*
     *  MARK: One-index transformations
     */

    /**
     *  Apply a one-index transformation and return the result.
     * 
     *  @param T            The basis transformation.
     * 
     *  @return The one-index-transformed one-electron operator.
     */
    DerivedOperator oneIndexTransformed(const Transformation& T) const {

        // Calculate the basis transformation for every component of the operator.
        const auto& parameters = this->allParameters();
        auto result = this->allParameters();

        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            const auto& f_i = parameters[i];

            result[i] = T.matrix().adjoint() * f_i + f_i * T.matrix();
        }

        return DerivedOperator {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};
    }


    /*
     *  MARK: Mulliken partitioning
     */

    /**
     *  Partition this one-electron operator according to the supplied Mulliken partitioning scheme.
     * 
     *  @param mulliken_partitioning                An encapsulation of the Mulliken partitioning scheme.
     * 
     *  @return A one-electron operator whose integrals/parameters/matrix elements correspond to the Mulliken-partitioning of this one-electron operator.
     */
    DerivedOperator partitioned(const MullikenPartitioning& mulliken_partitioning) const { return 0.5 * this->oneIndexTransformed(mulliken_partitioning.projectionMatrix()); }
};


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
struct OperatorTraits<SimpleSQOneElectronOperator<_Scalar, _Vectorizer, _DerivedOperator>> {

    // The scalar type used for a single parameter/matrix element/integral: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator that derives from `SimpleSQOneElectronOperator`, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = _DerivedOperator;
};


}  // namespace GQCP
