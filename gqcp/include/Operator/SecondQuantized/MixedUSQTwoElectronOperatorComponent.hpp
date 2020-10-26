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
#include "Mathematical/Representation/QCRankFourTensor.hpp"
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
    public SQOperatorStorage<QCRankFourTensor<_Scalar>, _Vectorizer, MixedUSQTwoElectronOperatorComponent<_Scalar, _Vectorizer>> {
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
    using SQOperatorStorage<QCRankFourTensor<_Scalar>, _Vectorizer, MixedUSQTwoElectronOperatorComponent<_Scalar, _Vectorizer>>::SQOperatorStorage;


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
    Self transformed(const UTransformationMatrixComponent<Scalar>& transformation_matrix, const Spin sigma) const override {

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

        return Self {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};  // TODO: Try to rewrite this.
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
 *  @tparam Scalar          The scalar type used for a single parameter: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
 */
template <typename Scalar, typename Vectorizer>
struct OperatorTraits<MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer>> {

    // A type that corresponds to the scalar version of the associated component of an unrestricted two-electron operator type.
    using ScalarOperator = ScalarMixedUSQTwoElectronOperatorComponent<Scalar>;

    // The type of transformation matrix that is naturally associated to a component of an unestricted two-electron operator.
    using TM = UTransformationMatrixComponent<Scalar>;

    // The type of the one-particle density matrix that is naturally associated a component of an unestricted two-electron operator.
    using OneDM = SpinResolved1DMComponent<Scalar>;

    // The type of the two-particle density matrix that is naturally associated a component of an unestricted two-electron operator.
    using TwoDM = SpinResolved2DMComponent<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar, typename Vectorizer>
struct BasisTransformableTraits<MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer>> {

    // The type of transformation matrix that is naturally associated to a component of an unestricted two-electron operator.
    using TM = UTransformationMatrixComponent<Scalar>;
};


}  // namespace GQCP
