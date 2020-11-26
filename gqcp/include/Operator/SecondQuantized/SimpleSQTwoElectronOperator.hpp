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
#include "Mathematical/Representation/SquareRankFourTensor.hpp"
#include "Mathematical/Representation/StorageArray.hpp"
#include "Operator/SecondQuantized/SQOperatorStorage.hpp"

#include <string>


namespace GQCP {


/**
 *  A second-quantized one-electron operator whose parameters are described by a single tensor.
 * 
 *  This class is used as a base class for `RSQTwoElectronOperator` and `GSQTwoElectronOperator`, since they both admit parameter representations using a single tensor, as opposed to `USQTwoElectronOperator`, which uses separate alpha- and beta- tensors. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                 The scalar type used for a single parameter/matrix element: real or complex.
 *  @tparam _Vectorizer             The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
 *  @tparam _DerivedOperator        The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
class SimpleSQTwoElectronOperator:
    public SQOperatorStorage<SquareRankFourTensor<_Scalar>, _Vectorizer, SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, _DerivedOperator>>,
    public BasisTransformable<_DerivedOperator>,
    public JacobiRotatable<_DerivedOperator> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = _DerivedOperator;

    // The type of 'this'.
    using Self = SimpleSQTwoElectronOperator<Scalar, Vectorizer, DerivedOperator>;

    // The matrix representation of the parameters of (one of the components of) the two-electron operator.
    using MatrixRepresentation = SquareRankFourTensor<Scalar>;

    // The type of one-electron operator that is naturally related to the derived two-electron operator.
    using DerivedSQOneElectronOperator = typename OperatorTraits<DerivedOperator>::SQOneElectronOperator;

    // The type of the one-particle density matrix that is naturally associated to the derived two-electron operator.
    using Derived1DM = typename OperatorTraits<DerivedOperator>::OneDM;

    // The type of the one-particle density matrix that is naturally associated to the derived two-electron operator.
    using Derived2DM = typename OperatorTraits<DerivedOperator>::TwoDM;

    // The type of transformation that is naturally associated to the derived two-electron operator.
    using Transformation = typename OperatorTraits<DerivedOperator>::Transformation;


private:
    // If this two-electron operator contains matrix elements that are modified to obey antisymmetry w.r.t. creation and annihilation indices.
    bool is_antisymmetrized = false;

    // If this two-electron operator contains matrix elements that are expressed as g_PQRS or (PQ|RS).
    bool is_expressed_using_chemists_notation = true;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SQOperatorStorage`'s constructors.
    using SQOperatorStorage<SquareRankFourTensor<Scalar>, Vectorizer, SimpleSQTwoElectronOperator<Scalar, Vectorizer, DerivedOperator>>::SQOperatorStorage;


    /*
     *  MARK: Calculations
     */

    /**
     *  Calculate the expectation value of this two-electron operator, given a two-electron density matrix. (This includes the prefactor 1/2.)
     * 
     *  @param d            The 2-DM (that represents the wave function).
     *
     *  @return The expectation values of all the components of the two-electron operator, with the given 2-DM.
     */
    StorageArray<Scalar, Vectorizer> calculateExpectationValue(const Derived2DM& d) const {

        if (this->numberOfOrbitals() != d.numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQTwoElectronOperator::calculateExpectationValue(const Derived2DM&): The given 2-DM's dimension is not compatible with the two-electron operator.");
        }


        // Calculate the expectation value for every component of the operator.
        const auto& parameters = this->allParameters();
        std::vector<Scalar> expectation_values(this->numberOfComponents());  // Zero-initialize the vector with a number of elements.

        for (size_t i = 0; i < this->numberOfComponents(); i++) {

            // Perform the actual contraction (with prefactor 0.5):
            //      0.5 g(p q r s) d(p q r s)
            Eigen::Tensor<Scalar, 0> contraction = 0.5 * parameters[i].template einsum<4>("pqrs, pqrs->", d);

            // As the contraction is a scalar (a tensor of rank 0), we should access using `operator(0)`.
            expectation_values[i] = contraction(0);
        }

        return StorageArray<Scalar, Vectorizer> {expectation_values, this->array.vectorizer()};  // convert std::array to Vector
    }


    /**
     *  Calculate the Fockian matrix for (each of the components of) this two-electron operator.
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
            throw std::invalid_argument("SimpleSQTwoElectronOperator::calculateFockianMatrix(const Derived1DM&, const Derived2DM&): The 1-DM's dimensions are not compatible with this two-electron operator.");
        }

        if (d.numberOfOrbitals() != this->numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQTwoElectronOperator::calculateFockianMatrix(const Derived1DM&, const Derived2DM&): The 2-DM's dimensions are not compatible with this two-electron operator.");
        }

        const auto& parameters = this->allParameters();                           // The parameters of this two-electron operator, as a vector.
        std::vector<SquareMatrix<double>> F_vector {this->numberOfComponents()};  // The resulting vector of Fockian matrices.

        // A KISS implementation of the calculation of the Fockian matrix.
        for (size_t i = 0; i < this->numberOfComponents(); i++) {

            const auto& g_i = parameters[i];  // The matrix representation of the parameters of the i-th component.

            // Calculate the Fockian matrix for every component and add it to the vector.
            SquareMatrix<Scalar> F_i = SquareMatrix<double>::Zero(this->numberOfOrbitals());
            for (size_t p = 0; p < this->numberOfOrbitals(); p++) {
                for (size_t q = 0; q < this->numberOfOrbitals(); q++) {

                    for (size_t r = 0; r < this->numberOfOrbitals(); r++) {
                        for (size_t s = 0; s < this->numberOfOrbitals(); s++) {
                            for (size_t t = 0; t < this->numberOfOrbitals(); t++) {
                                F_i(p, q) += 0.5 * g_i(q, r, s, t) * (d(p, r, s, t) + d(r, p, s, t));  // Include a factor 1/2 to accommodate for response density matrices.
                            }
                        }
                    }
                }
            }  // F_i elements loop
            F_vector[i] = F_i;
        }

        return StorageArray<SquareMatrix<double>, Vectorizer> {F_vector, this->array.vectorizer()};
    }


    /**
     *  Calculate the super-Fockian matrix for (each of the components of) this two-electron operator.
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
            throw std::invalid_argument("SimpleSQTwoElectronOperator::calculateFockianMatrix(const Derived1DM&, const Derived2DM&): The given 1-DM's dimensions are not compatible with this two-electron operator.");
        }

        if (d.numberOfOrbitals() != this->numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQTwoElectronOperator::calculateFockianMatrix(const Derived1DM&, const Derived2DM&): The given 2-DM's dimensions are not compatible with this two-electron operator.");
        }


        const auto& parameters = this->allParameters();                                   // The parameters of this two-electron operator, as a vector.
        std::vector<SquareRankFourTensor<double>> G_vector {this->numberOfComponents()};  // The resulting vector of super-Fockian matrices.
        const auto F_vector = this->calculateFockianMatrix(D, d).elements();              // The Fockian matrices are necessary in the calculation of the super-Fockian matrices.


        // A KISS implementation of the calculation of the super-Fockian matrix
        for (size_t i = 0; i < this->numberOfComponents(); i++) {

            const auto& g_i = parameters[i];  // The matrix representation of the parameters of the i-th component.
            const auto& F_i = F_vector[i];    // The Fockian matrix of the i-th component.

            // Calculate the super-Fockian matrix for every component and add it to the array. Add factors 1/2 to accommodate for response density matrices.
            SquareRankFourTensor<Scalar> G_i {this->numberOfOrbitals()};
            G_i.setZero();
            for (size_t p = 0; p < this->numberOfOrbitals(); p++) {
                for (size_t q = 0; q < this->numberOfOrbitals(); q++) {
                    for (size_t r = 0; r < this->numberOfOrbitals(); r++) {
                        for (size_t s = 0; s < this->numberOfOrbitals(); s++) {

                            if (q == r) {
                                G_i(p, q, r, s) += F_i(p, s);
                            }

                            for (size_t t = 0; t < this->numberOfOrbitals(); t++) {
                                for (size_t u = 0; u < this->numberOfOrbitals(); u++) {
                                    double value = g_i(s, t, q, u) * (d(r, t, p, u) + d(t, r, u, p));
                                    value -= g_i(s, t, u, p) * (d(r, t, u, q) + d(t, r, q, u));
                                    value -= g_i(s, p, t, u) * (d(r, q, t, u) + d(q, r, u, t));

                                    G_i(p, q, r, s) += 0.5 * value;
                                }
                            }
                        }
                    }
                }
            }  // G_i elements loop
            G_vector[i] = G_i;
        }

        return StorageArray<SquareRankFourTensor<double>, Vectorizer> {G_vector, this->array.vectorizer()};
    }


    /**
     *  @return The one-electron operator that is the difference between this two-electron operator (E_PQRS) and a product of one-electron operators (E_PQ E_RS).
     */
    DerivedSQOneElectronOperator effectiveOneElectronPartition() const {

        // Prepare some variables.
        const auto& g = this->allParameters();  // The parameters of this two-electron operator, as a vector.

        const auto K = this->numberOfOrbitals();
        DerivedSQOneElectronOperator k_op = DerivedSQOneElectronOperator::Zero(K);  // 'op' for operator.
        auto& k = k_op.allParameters();

        // A KISS-implementation of the elements of the one-electron partition.
        for (size_t i = 0; i < this->numberOfComponents(); i++) {

            for (size_t p = 0; p < K; p++) {
                for (size_t q = 0; q < K; q++) {
                    for (size_t r = 0; r < K; r++) {
                        k[i](p, q) -= 0.5 * g[i](p, r, r, q);
                    }
                }
            }
        }

        return k_op;
    }


    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the resulting two-electron integrals.
     * 
     *  @param T            The basis transformation
     * 
     *  @return The basis-transformed one-electron integrals.
     */
    DerivedOperator transformed(const Transformation& T) const override {

        // Since we're only getting T as a matrix, we should convert it to an appropriate tensor to perform contractions.
        // Although not a necessity for the einsum implementation, it makes it a lot easier to follow the formulas.
        const GQCP::Tensor<Scalar, 2> T_tensor = Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>>(T.matrix().data(), T.matrix().rows(), T.matrix().cols());

        // We calculate the conjugate as a tensor as well.
        const GQCP::Tensor<Scalar, 2> T_conjugate = T_tensor.conjugate();

        // Calculate the basis transformation for every component of the operator.
        const auto& parameters = this->allParameters();
        auto result = this->allParameters();

        for (size_t i = 0; i < this->numberOfComponents(); i++) {

            // We will have to do four single contractions
            // g(T U V W)  T^*(V R) -> a(T U R W)
            // a(T U R W)  T(W S) -> b(T U R S)
            const auto temp_1 = parameters[i].template einsum<1>("TUVW,VR->TURW", T_conjugate).template einsum<1>("TURW,WS->TURS", T_tensor);

            // T(U Q)  b(T U R S) -> c(T Q R S)
            const auto temp_2 = T_tensor.template einsum<1>("UQ,TURS->TQRS", temp_1);

            // T^*(T P)  c(T Q R S) -> g'(P Q R S)
            const auto transformed_component = T_conjugate.template einsum<1>("TP,TQRS->PQRS", temp_2);

            result[i] = transformed_component;
        }

        return DerivedOperator {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};
    }


    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<DerivedOperator>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedOperator>::rotated;


    /*
     *  MARK: Conforming to `JacobiRotatable`
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_rotation          The Jacobi rotation.
     * 
     *  @return The Jacobi-transformed object.
     */
    DerivedOperator rotated(const JacobiRotation& jacobi_rotation) const override {

        // While waiting for an analogous Eigen::Tensor Jacobi module, we implement this rotation by constructing a Jacobi rotation matrix and then simply doing a rotation with it.
        const auto J = Transformation::FromJacobi(jacobi_rotation, this->numberOfOrbitals());
        return this->rotated(J);
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<DerivedOperator>::rotate;


    /*
     *  MARK: Antisymmetrizing
     */

    /**
     *  @return If these two-electron integrals are considered to be antisymmetrized.
     * 
     *  @note If so, these integrals represent:
     *      - if they are expressed using chemist's notation:       g_{PQRS} - g_{PSRQ}, i.e. they are antisymmetric upon interchanging the indices PR or QS
     *      - if they are expressed using physicist's notation:     <PQ|RS> - <PQ|SR>, i.e. they are antisymmetric upon interchanging the indices PQ or RS
     */
    bool isAntisymmetrized() const { return this->is_antisymmetrized; }


    /**
     *  @return An antisymmetrized version of this two-electron operator, i.e. one with matrix elements that are modified to obey antisymmetry w.r.t. creation and annihilation indices.
     * 
     *  @note If the integrals are expressed using
     *      - chemist's notation g_{PQRS}, return g_{PQRS} - g_{PSRQ}
     *      - physicist's notation <PQ|RS>, return <PQ||RS> = <PQ|RS> - <PQ|SR>
     */
    Self antisymmetrized() const {

        // Attempt to modify a copy of these integrals if they haven't been antisymmetrized already.
        const auto& parameters = this->allParameters();
        auto copy = *this;
        auto& copy_parameters = copy.allParameters();

        if (!(this->isAntisymmetrized())) {

            // Determine the shuffle-indices for the matrix elements that will be subtracted.
            Eigen::array<int, 4> shuffle_indices;
            if (this->isExpressedUsingChemistsNotation()) {
                shuffle_indices = Eigen::array<int, 4> {0, 3, 2, 1};
            } else {  // Expressed using physicist's notation.
                shuffle_indices = Eigen::array<int, 4> {0, 1, 3, 2};
            }

            auto& copy_parameters = copy.allParameters();
            for (size_t i = 0; i < this->numberOfComponents(); i++) {
                copy_parameters[i] -= parameters[i].shuffle(shuffle_indices);
            }

            copy.is_antisymmetrized = true;
        }

        return copy;
    }


    /**
     *  In-place antisymmetrize this two-electron operator.
     * 
     *  @note If the integrals are expressed using
     *      - chemist's notation g_{PQRS}, they are modified to g_{PQRS} - g_{PSRQ}
     *      - physicist's notation <PQ|RS>, they are modified to <PQ||RS> = <PQ|RS> - <PQ|SR>
     */
    void antisymmetrize() { *this = this->antisymmetrized(); }


    /*
     *  MARK: Chemist's or physicist's notation
     */

    /**
     *  @return If this two-electron operator's integrals are expressed using chemist's notation g_{PQRS}, i.e. (PQ|RS)
     */
    bool isExpressedUsingChemistsNotation() const { return this->is_expressed_using_chemists_notation; }

    /**
     *  @return If this two-electron operator's integrals are expressed using physicist's notation <PQ|RS>
     */
    bool isExpressedUsingPhysicistsNotation() const { return !(this->isExpressedUsingChemistsNotation()); }

    /**
     *  @return The two-electron operator with integrals changed to chemist's notation (from physicist's notation).
     */
    Self convertedToChemistsNotation() const {

        // Attempt to modify a copy if these integrals are expressed in physicist's notation.
        const auto& parameters = this->allParameters();
        auto copy = *this;
        auto& copy_parameters = copy.allParameters();

        if (this->isExpressedUsingPhysicistsNotation()) {

            Eigen::array<int, 4> shuffle_indices {0, 2, 1, 3};

            for (size_t i = 0; i < this->numberOfComponents(); i++) {
                copy_parameters[i] = SquareRankFourTensor<double>(parameters[i].shuffle(shuffle_indices));
            }

            copy.is_expressed_using_chemists_notation = true;
        }
        return copy;
    }


    /**
     *  @return The two-electron operator with integrals changed to physicist's notation (from chemist's notation).
     */
    Self convertedToPhysicistsNotation() const {

        // Attempt to modify a copy if these integrals are expressed in chemist's notation.
        const auto& parameters = this->allParameters();
        auto copy = *this;
        auto& copy_parameters = copy.allParameters();

        if (this->isExpressedUsingChemistsNotation()) {

            Eigen::array<int, 4> shuffle_indices {0, 2, 1, 3};

            for (size_t i = 0; i < this->numberOfComponents(); i++) {
                copy_parameters[i] = SquareRankFourTensor<double>(parameters[i].shuffle(shuffle_indices));
            }

            copy.is_expressed_using_chemists_notation = false;
        }
        return copy;
    }


    /**
     *  In-place change this two-electron operator's integrals to chemist's notation (from physicist's notation).
     */
    void convertToChemistsNotation() { *this = this->convertedToChemistsNotation(); }

    /**
     *  In-place change this two-electron operator's integrals to physicist's notation (from chemist's notation).
     */
    void convertToPhysicistsNotation() { *this = this->convertedToPhysicistsNotation(); }
};


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
struct OperatorTraits<SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, _DerivedOperator>> {

    // The scalar type used for a single parameter/matrix element/integral: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
    using Vectorizer = _Vectorizer;

    // The type of the operator that derives from `SimpleSQTwoElectronOperator`, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = _DerivedOperator;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
struct BasisTransformableTraits<SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, _DerivedOperator>> {

    // The type of the transformation for which the basis transformation should be defined.
    using Transformation = typename OperatorTraits<_DerivedOperator>::Transformation;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
struct JacobiRotatableTraits<SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, _DerivedOperator>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
