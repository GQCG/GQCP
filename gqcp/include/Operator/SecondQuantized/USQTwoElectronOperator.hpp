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


#include "Basis/Transformations/DoublySpinResolvedBasisTransformable.hpp"
#include "Basis/Transformations/DoublySpinResolvedJacobiRotatable.hpp"
#include "Basis/Transformations/UJacobiRotation.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "DensityMatrix/SpinResolved1DM.hpp"
#include "DensityMatrix/SpinResolved2DM.hpp"
#include "Mathematical/Functions/VectorSpaceArithmetic.hpp"
#include "Operator/SecondQuantized/MixedUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/PureUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "QuantumChemical/DoublySpinResolvedBase.hpp"
#include "QuantumChemical/spinor_tags.hpp"


namespace GQCP {


/**
 *  A class that represents an 'unrestricted second-quantized two-electron operator' suitable for the projection of the non-relativistic Hamiltonian onto an unrestricted spinor basis. It holds the tensor representation of its parameters for both spin components and both mixed spin components, which are (usually) integrals over first-quantized operators.
 * 
 *  @tparam _Scalar             The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer         The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
 */
template <typename _Scalar, typename _Vectorizer>
class USQTwoElectronOperator:
    public DoublySpinResolvedBase<PureUSQTwoElectronOperatorComponent<_Scalar, _Vectorizer>, MixedUSQTwoElectronOperatorComponent<_Scalar, _Vectorizer>, USQTwoElectronOperator<_Scalar, _Vectorizer>>,
    public DoublySpinResolvedBasisTransformable<USQTwoElectronOperator<_Scalar, _Vectorizer>>,
    public DoublySpinResolvedJacobiRotatable<USQTwoElectronOperator<_Scalar, _Vectorizer>>,
    public VectorSpaceArithmetic<USQTwoElectronOperator<_Scalar, _Vectorizer>, _Scalar> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
    using Vectorizer = _Vectorizer;

    // The type of 'this'.
    using Self = USQTwoElectronOperator<Scalar, Vectorizer>;

    // The type of the transformation that is naturally related to an unrestricted two-electron operator.
    using Transformation = UTransformation<Scalar>;

    // The type of Jacobi rotation that is naturally related to an unrestricted two-electron operator.
    using JacobiRotationType = UJacobiRotation;

    // The spinor tag corresponding to a `USQTwoElectronOperator`.
    using SpinorTag = UnrestrictedSpinOrbitalTag;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `DoublySpinResolvedBase`'s constructors.
    using DoublySpinResolvedBase<PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer>, MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer>, USQTwoElectronOperator<Scalar, Vectorizer>>::DoublySpinResolvedBase;


    /**
     *  Create an `USQTwoElectronOperator` from all the matrix representations of its components.
     * 
     *  @param gs_aa                All the matrix representations of the components of the alpha-alpha-part of the unrestricted two-electron operator.
     *  @param gs_ab                All the matrix representations of the components of the alpha-beta-part of the unrestricted two-electron operator.
     *  @param gs_ba                All the matrix representations of the components of the beta-alpha-part of the unrestricted two-electron operator.
     *  @param gs_bb                All the matrix representations of the components of the beta-beta-part of the unrestricted two-electron operator.
     *  @param vectorizer           The type of the vectorizer that relates a one-dimensional storage of matrix representations to the tensor structure of the second-quantized operator. Equal for all the different spin-components.
     * 
     *  @tparam N                   The number of components for the different spin components of the unrestricted two-electron operator.
     */
    template <size_t N>
    USQTwoElectronOperator(const std::array<SquareRankFourTensor<Scalar>, N>& gs_aa, const std::array<SquareRankFourTensor<Scalar>, N>& gs_ab, const std::array<SquareRankFourTensor<Scalar>, N>& gs_ba, const std::array<SquareRankFourTensor<Scalar>, N>& gs_bb, const Vectorizer& vectorizer) :

        // Encapsulate the array of matrix representations in the different spin-components components, and put them together to form the `USQTwoElectronOperator`.
        DoublySpinResolvedBase<PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer>, MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer>, USQTwoElectronOperator<Scalar, Vectorizer>>(
            PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer> {StorageArray<SquareRankFourTensor<Scalar>, Vectorizer> {gs_aa, vectorizer}},
            MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer> {StorageArray<SquareRankFourTensor<Scalar>, Vectorizer> {gs_ab, vectorizer}},
            MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer> {StorageArray<SquareRankFourTensor<Scalar>, Vectorizer> {gs_ba, vectorizer}},
            PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer> {StorageArray<SquareRankFourTensor<Scalar>, Vectorizer> {gs_bb, vectorizer}}) {

        // Check if the given tensor representations have the same dimensions, for each spin part.
        const auto dimension_of_first_aa = gs_aa[0].dimension();
        const auto dimension_of_first_ab = gs_ab[0].dimension();
        const auto dimension_of_first_ba = gs_ba[0].dimension();
        const auto dimension_of_first_bb = gs_bb[0].dimension();

        for (size_t i = 1; i < N; i++) {

            const auto dimension_of_ith_aa = gs_aa[i].dimension();
            const auto dimension_of_ith_ab = gs_ab[i].dimension();
            const auto dimension_of_ith_ba = gs_ba[i].dimension();
            const auto dimension_of_ith_bb = gs_bb[i].dimension();

            if ((dimension_of_first_aa != dimension_of_ith_aa) || (dimension_of_first_ab != dimension_of_ith_ab) || (dimension_of_first_ba != dimension_of_ith_ba) || (dimension_of_first_bb != dimension_of_ith_bb)) {
                throw std::invalid_argument("USQTwoElectronOperator(const std::array<SquareRankFourTensor<Scalar>, N>&, const std::array<SquareRankFourTensor<Scalar>, N>&, const std::array<SquareRankFourTensor<Scalar>, N>&, const std::array<SquareRankFourTensor<Scalar>, N>&): The given tensor representations do not have the same dimensions for either the alpha, beta or one of the mixed components.");
            }
        }
    }


    /**
     *  A constructor for `ScalarUSQTwoElectronOperator`s that doesn't require the argument to be an array of just one element.
     * 
     *  @param g_aa           The tensor representation of the alpha-alpha integrals of the scalar second-quantized operator.
     *  @param g_ab           The tensor representation of the alpha-beta integrals of the scalar second-quantized operator.
     *  @param g_ba           The tensor representation of the beta-alpha integrals of the scalar second-quantized operator.
     *  @param g_bb           The tensor representation of the beta-beta integrals of the scalar second-quantized operator.
     * 
     *  @note This constructor is only available for `ScalarUSQTwoElectronOperators` (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <typename Z = Vectorizer>
    USQTwoElectronOperator(const SquareRankFourTensor<Scalar>& g_aa, const SquareRankFourTensor<Scalar>& g_ab, const SquareRankFourTensor<Scalar>& g_ba, const SquareRankFourTensor<Scalar>& g_bb,
                           typename std::enable_if<std::is_same<Z, ScalarVectorizer>::value>::type* = 0) :
        USQTwoElectronOperator(std::array<SquareRankFourTensor<Scalar>, 1> {g_aa}, std::array<SquareRankFourTensor<Scalar>, 1> {g_ab}, std::array<SquareRankFourTensor<Scalar>, 1> {g_ba}, std::array<SquareRankFourTensor<Scalar>, 1> {g_bb}, ScalarVectorizer()) {}


    /**
     *  The default constructor.
     */
    USQTwoElectronOperator() :
        USQTwoElectronOperator(USQTwoElectronOperator<Scalar, Vectorizer>::Zero(0)) {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Construct an unrestricted two-electron operator with parameters that are zero. The dimensions different spin components are equal.
     * 
     *  @param dim        The dimension of the axes of the the matrix representation of the different spin components. Equal for all spin components.
     */
    static USQTwoElectronOperator<Scalar, Vectorizer> Zero(const size_t dim) {

        const auto zero_component_pure = PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer>::Zero(dim);
        const auto zero_component_mixed = MixedUSQTwoElectronOperatorComponent<Scalar, Vectorizer>::Zero(dim);

        return USQTwoElectronOperator<Scalar, Vectorizer> {zero_component_pure, zero_component_mixed, zero_component_mixed, zero_component_pure};
    }


    /**
     *  Construct an `USQTwoElectronOperator` from an `RSQTwoElectronOperator.
     * 
     *  @param g_restricted             The restricted two-electron operator that should be converted.
     */
    static USQTwoElectronOperator<Scalar, Vectorizer> FromRestricted(const RSQTwoElectronOperator<Scalar, Vectorizer>& g_restricted) {

        return USQTwoElectronOperator<Scalar, Vectorizer> {g_restricted.alphaAlpha(), g_restricted.alphaBeta(), g_restricted.betaAlpha(), g_restricted.betaBeta()};
    }


    /*
     *  MARK: Calculations
     */

    /**
     *  Calculate the expectation value of this two-electron operator.
     * 
     *  @param d                    The spin-resolved 2-DM (that represents the wave function).
     *
     *  @return the expectation values of all the components of the two-electron operator, with the given 2-DMs: this includes the prefactor 1/2
     */
    StorageArray<Scalar, Vectorizer> calculateExpectationValue(const SpinResolved2DM<Scalar>& d) const {

        // Calculate the sum of all spin contributions.
        // Unfortunately, we can't give `StorageArray` out-of-the-box vector-space arithmetic (since the default scalar type for scalar multiplication is unknown), so we'll have to do the summation ourselves.
        const auto expectation_value_elements_aa = this->alphaAlpha().calculateExpectationValue(d.alphaAlpha()).elements();
        const auto expectation_value_elements_ab = this->alphaBeta().calculateExpectationValue(d.alphaBeta()).elements();
        const auto expectation_value_elements_ba = this->betaAlpha().calculateExpectationValue(d.betaAlpha()).elements();
        const auto expectation_value_elements_bb = this->betaBeta().calculateExpectationValue(d.betaBeta()).elements();

        auto result_elements = expectation_value_elements_aa;
        const auto vectorizer = this->alphaAlpha().vectorizer();
        for (size_t i = 0; i < vectorizer.numberOfElements(); i++) {
            result_elements[i] += expectation_value_elements_ab[i] + expectation_value_elements_ba[i] + expectation_value_elements_bb[i];
        }

        return StorageArray<Scalar, Vectorizer> {result_elements, vectorizer};
    }


    /*
     *  MARK: General information
     */

    /**
     *  @param sigma            Alpha or beta.
     *  @param tau              Alpha or beta.
     * 
     *  @return The number of orbitals (spinors or spin-orbitals, depending on the context) that are related to the sigma-tau part of the unrestricted two-electron operator.
     */
    size_t numberOfOrbitals(const Spin sigma, const Spin tau) const {

        if (sigma == Spin::alpha && tau == Spin::beta) {
            return this->alphaAlpha().numberOfOrbitals();
        } else if (sigma == Spin::alpha && tau == Spin::beta) {
            return this->alphaBeta().numberOfOrbitals();
        } else if (sigma == Spin::beta && tau == Spin::alpha) {
            return this->betaAlpha().numberOfOrbitals();
        } else {
            return this->betaBeta().numberOfOrbitals();
        }
    }


    /*
     *  @return The number of orbital related to the alpha-alpha part of the unrestricted two-electron operator.
     * 
     *  @note It is advised to only use this API when it is known that all spin-components of the two-electron operator are equal.
     */
    size_t numberOfOrbitals() const { return this->alphaAlpha().numberOfOrbitals(); }


    /**
     *  MARK: Enabling basis transformations
     */

    // Since `rotate` and `rotated` are both defined in `DoublySpinResolvedBasisTransformable` and `DoublySpinResolvedJacobiRotatable`, we have to explicitly enable these methods here.

    // Allow the `rotate` method from `DoublySpinResolvedBasisTransformable`, since there's also a `rotate` from `DoublySpinResolvedJacobiRotatable`.
    using DoublySpinResolvedBasisTransformable<Self>::rotate;

    // Allow the `rotated` method from `DoublySpinResolvedBasisTransformable`, since there's also a `rotated` from `DoublySpinResolvedJacobiRotatable`.
    using DoublySpinResolvedBasisTransformable<Self>::rotated;

    // Allow the `rotate` method from `DoublySpinResolvedJacobiRotatable`, since there's also a `rotate` from `DoublySpinResolvedBasisTransformable`.
    using DoublySpinResolvedJacobiRotatable<Self>::rotate;

    // Allow the `rotated` method from `DoublySpinResolvedJacobiRotatable`, since there's also a `rotated` from `DoublySpinResolvedBasisTransformable`.
    using DoublySpinResolvedJacobiRotatable<Self>::rotated;


    /*
     *  MARK: Conforming to VectorSpaceArithmetic
     */

    /**
     *  Addition-assignment.
     */
    Self& operator+=(const Self& rhs) override {

        // Add the spin-components.
        this->alphaAlpha() += rhs.alphaAlpha();
        this->alphaBeta() += rhs.alphaBeta();
        this->betaAlpha() += rhs.betaAlpha();
        this->betaBeta() += rhs.betaBeta();

        return *this;
    }


    /**
     *  Scalar multiplication-assignment, to be implemented in derived classes.
     */
    Self& operator*=(const Scalar& a) override {

        // Multiply the spin-components with the scalar.
        this->alphaAlpha() *= a;
        this->alphaBeta() *= a;
        this->betaAlpha() *= a;
        this->betaBeta() *= a;

        return *this;
    }
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like USQTwoElectronOperator, i.e. with scalar-like access.
template <typename Scalar>
using ScalarUSQTwoElectronOperator = USQTwoElectronOperator<Scalar, ScalarVectorizer>;

// A vector-like USQTwoElectronOperator, i.e. with vector-like access.
template <typename Scalar>
using VectorUSQTwoElectronOperator = USQTwoElectronOperator<Scalar, VectorVectorizer>;

// A matrix-like USQTwoElectronOperator, i.e. with matrix-like access.
template <typename Scalar>
using MatrixUSQTwoElectronOperator = USQTwoElectronOperator<Scalar, MatrixVectorizer>;

// A tensor-like USQTwoElectronOperator, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorUSQTwoElectronOperator = USQTwoElectronOperator<Scalar, TensorVectorizer<N>>;


/**
 *  A type that provides compile-time information (traits) on `USQTwoElectronOperator` that is otherwise not accessible through a public class alias.
 * 
 *  @tparam Scalar          The scalar type used for a single parameter: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of two-electron operators. This distinction is carried over from `USQOneElectronOperator`.
 */
template <typename Scalar, typename Vectorizer>
struct OperatorTraits<USQTwoElectronOperator<Scalar, Vectorizer>> {

    // A type that corresponds to the scalar version of the associated unrestricted one-electron operator type.
    using ScalarOperator = ScalarUSQTwoElectronOperator<Scalar>;

    // The type of transformation that is naturally associated to an unrestricted two-electron operator.
    using Transformation = UTransformation<Scalar>;

    // The type of the one-particle density matrix that is naturally associated an unrestricted two-electron operator.
    using OneDM = SpinResolved1DM<Scalar>;

    // The type of the two-particle density matrix that is naturally associated an unrestricted two-electron operator.
    using TwoDM = SpinResolved2DM<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar, typename Vectorizer>
struct BasisTransformableTraits<USQTwoElectronOperator<Scalar, Vectorizer>> {

    // The type of transformation that is naturally related to a `USQTwoElectronOperator`.
    using Transformation = UTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar, typename Vectorizer>
struct JacobiRotatableTraits<USQTwoElectronOperator<Scalar, Vectorizer>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = UJacobiRotation;
};


}  // namespace GQCP
