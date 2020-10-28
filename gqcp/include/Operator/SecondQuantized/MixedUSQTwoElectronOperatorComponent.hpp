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


#include "Basis/Transformations/UTransformationMatrixComponent.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"
#include "DensityMatrix/SpinResolved2DMComponent.hpp"
#include "Mathematical/Representation/SquareRankFourTensor.hpp"
#include "Operator/SecondQuantized/SQOperatorStorage.hpp"
#include "QuantumChemical/Spin.hpp"


namespace GQCP {


/**
 *  One of the pure (i.e. alpha-alpha or beta-beta) spin components of an unrestricted two-electron operator.
 * 
 *  @tparam _Scalar         The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
 */
template <typename _Scalar, typename _Vectorizer>
class MixedUSQTwoElectronOperatorComponent:
    public SQOperatorStorage<SquareRankFourTensor<_Scalar>, _Vectorizer, MixedUSQTwoElectronOperatorComponent<_Scalar, _Vectorizer>> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
    using Vectorizer = _Vectorizer;

    // The type of 'this'.
    using Self = MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SQOperatorStorage`'s constructors.
    using SQOperatorStorage<SquareRankFourTensor<_Scalar>, _Vectorizer, MixedUSQTwoElectronOperatorComponent<_Scalar, _Vectorizer>>::SQOperatorStorage;


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
    StorageArray<Scalar, Vectorizer> calculateExpectationValue(const SpinResolved2DMComponent<Scalar>& d) const {

        // FIXME: This is duplicate code from `SimpleSQTwoElectronOperator`. It would be double work to introduce a temporary intermediate class before resolving issue #559 (https://github.com/GQCG/GQCP/issues/559), which is why we chose to just copy-paste the implementation.

        if (this->numberOfOrbitals() != d.numberOfOrbitals()) {
            throw std::invalid_argument("MixedUSQTwoElectronOperatorComponent::calculateExpectationValue(const SpinResolved2DMComponent&): The given 2-DM's dimension is not compatible with the two-electron operator.");
        }


        // Calculate the expectation value for every component of the operator.
        const auto& parameters = this->allParameters();
        std::vector<Scalar> expectation_values(this->numberOfComponents());  // Zero-initialize the vector with a number of elements.

        for (size_t i = 0; i < this->numberOfComponents(); i++) {

            // Specify the contractions for the relevant contraction of the two-electron integrals/parameters/matrix elements and the 2-DM:
            //      0.5 g(p q r s) d(p q r s)
            Eigen::array<Eigen::IndexPair<int>, 4> contractions {Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(1, 1), Eigen::IndexPair<int>(2, 2), Eigen::IndexPair<int>(3, 3)};

            // Perform the actual contraction.
            Eigen::Tensor<Scalar, 0> contraction = 0.5 * parameters[i].contract(d.Eigen(), contractions);

            // As the contraction is a scalar (a tensor of rank 0), we should access using `operator(0)`.
            expectation_values[i] = contraction(0);
        }

        return StorageArray<Scalar, Vectorizer> {expectation_values, this->array.vectorizer()};  // convert std::array to Vector
    }


    /*
     *  MARK: basis transformations
     */

    /**
     *  Apply the basis transformation for the spin component sigma, and return the resulting two-electron integrals.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     *  @param sigma                        Alpha indicates a transformation of the first two axes, beta indicates a transformation of the second two axes.
     * 
     *  @return The basis-transformed two-electron operator.
     * 
     *  @note We apologize for this half-baked API. It is currently present in the code, while issue #559 (https://github.com/GQCG/GQCP/issues/688) is being implemented.
     */
    Self transformed(const UTransformationMatrixComponent<Scalar>& transformation_matrix, const Spin sigma) const {

        // // Since we're only getting T as a matrix, we should convert it to an appropriate tensor to perform contractions.
        // const Tensor<double, 2> T = Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>>(transformation_matrix.data(), transformation_matrix.rows(), transformation_matrix.cols());

        // const Tensor<double, 2> T_conjugate = T.conjugate();


        // // Depending on the given spin-component, we should either transform the first two, or the second two axes.
        // auto result = *this;

        // for (size_t i = 0; i < this->numberOfComponents(); i++) {
        //     switch (sigma) {
        //     case Spin::alpha: {
        //         result.allParameters(i)
        //             .template einsum<1>(T_conjugate, "PQRS", "PT", "TQRS")
        //             .template einsum<1>(T, "TQRS", "QU", "TURS");
        //     }

        //     case Spin::beta: {
        //         result.allParameters(i)
        //             .template einsum<1>(T_conjugate, "PQRS", "RT", "PQTS")
        //             .template einsum<1>(T, "PQTS", "SU", "PQTU");
        //     }
        //     }
        // }

        // return result;

        // Since we're only getting T as a matrix, we should convert it to an appropriate tensor to perform contractions.
        // const Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> T_tensor {transformation_matrix.data(), transformation_matrix.rows(), transformation_matrix.cols()};

        const size_t first_contraction_index = 2 * sigma;
        const size_t second_contraction_index = 2 * sigma + 1;

        auto result = this->allParameters();

        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            result[i].template contractWithMatrix<Scalar>(transformation_matrix, first_contraction_index);
            result[i].template contractWithMatrix<Scalar>(transformation_matrix, second_contraction_index);
        }

        return Self {StorageArray<SquareRankFourTensor<Scalar>, Vectorizer>(result, this->array.vectorizer())};  // TODO: Try to rewrite this.
    }


    /**
     *  In-place apply the basis transformation for the spin component sigma.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     *  @param sigma                        Alpha indicates a transformation of the first two axes, beta indicates a transformation of the second two axes.
     * 
     *  @note We apologize for this half-baked API. It is currently present in the code, while issue #559 (https://github.com/GQCG/GQCP/issues/688) is being implemented.
     */
    void transform(const UTransformationMatrixComponent<Scalar>& transformation_matrix, const Spin sigma) {

        auto result = this->transformed(transformation_matrix, sigma);
        *this = result;
    }


    /**
     *  Apply the basis rotation for the spin component sigma, and return the resulting two-electron integrals.
     * 
     *  @param jacobi_rotation              The Jacobi rotation.
     *  @param sigma                        Alpha indicates a transformation of the first two axes, beta indicates a transformation of the second two axes.
     * 
     *  @return The basis-rotated two-electron operator.
     * 
     *  @note We apologize for this half-baked API. It is currently present in the code, while issue #559 (https://github.com/GQCG/GQCP/issues/688) is being implemented.
     */
    Self rotated(const JacobiRotation& jacobi_rotation, const Spin sigma) const {

        const auto J = UTransformationMatrixComponent<Scalar>::FromJacobi(jacobi_rotation, this->numberOfOrbitals());
        return this->transformed(J, sigma);
    }


    /**
     *  In-place apply the basis rotation for the spin component sigma, and return the resulting two-electron integrals.
     * 
     *  @param jacobi_rotation              The Jacobi rotation.
     *  @param sigma                        Alpha indicates a transformation of the first two axes, beta indicates a transformation of the second two axes.
     * 
     *  @note We apologize for this half-baked API. It is currently present in the code, while issue #559 (https://github.com/GQCG/GQCP/issues/688) is being implemented.
     */
    void rotate(const JacobiRotation& jacobi_rotation, const Spin sigma) {

        auto result = this->rotated(jacobi_rotation, sigma);
        *this = result;
    }
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like MixedUSQTwoElectronOperatorComponent, i.e. with scalar-like access.
template <typename Scalar>
using ScalarMixedUSQTwoElectronOperatorComponent = MixedUSQTwoElectronOperatorComponent<Scalar, ScalarVectorizer>;

// A vector-like MixedUSQTwoElectronOperatorComponent, i.e. with vector-like access.
template <typename Scalar>
using VectorMixedUSQTwoElectronOperatorComponent = MixedUSQTwoElectronOperatorComponent<Scalar, VectorVectorizer>;

// A matrix-like MixedUSQTwoElectronOperatorComponent, i.e. with matrix-like access.
template <typename Scalar>
using MatrixMixedUSQTwoElectronOperatorComponent = MixedUSQTwoElectronOperatorComponent<Scalar, MatrixVectorizer>;

// A tensor-like MixedUSQTwoElectronOperatorComponent, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorMixedUSQTwoElectronOperatorComponent = MixedUSQTwoElectronOperatorComponent<Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information (traits) on `MixedUSQTwoElectronOperatorComponent` that is otherwise not accessible through a public class alias.
 * 
 *  @tparam _Scalar         The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
 */
template <typename _Scalar, typename _Vectorizer>
struct OperatorTraits<MixedUSQTwoElectronOperatorComponent<_Scalar, _Vectorizer>> {

    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
    using Vectorizer = _Vectorizer;

    // A type that corresponds to the scalar version of the associated component of an unrestricted two-electron operator type.
    using ScalarOperator = ScalarMixedUSQTwoElectronOperatorComponent<Scalar>;

    // The type of the final derived operator (that derives from SQOperatorStorage), enabling CRTP and compile-time polymorphism. VectorSpaceArithmetic (and other functionality) should be implemented on the **final** deriving class, not on intermediate classes.
    using DerivedOperator = MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer>;

    // The type of transformation matrix that is naturally associated to a component of an unestricted two-electron operator.
    using TM = UTransformationMatrixComponent<Scalar>;

    // The type of the one-particle density matrix that is naturally associated a component of an unestricted two-electron operator.
    using OneDM = SpinResolved1DMComponent<Scalar>;

    // The type of the two-particle density matrix that is naturally associated a component of an unestricted two-electron operator.
    using TwoDM = SpinResolved2DMComponent<Scalar>;
};


}  // namespace GQCP
