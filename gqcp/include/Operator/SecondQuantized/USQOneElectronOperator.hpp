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


#include "Basis/Transformations/SpinResolvedBasisTransformable.hpp"
#include "Basis/Transformations/UTransformationMatrix.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"
#include "QuantumChemical/spinor_tags.hpp"


namespace GQCP {


/**
 *  A class that represents an 'unrestricted second-quantized one-electron operator'. This type of operator is suitable for the projection of the non-relativistic Hamiltonian onto an unrestricted spinor basis. It holds the matrix representation of its parameters for both spin components.
 *
 *  @tparam _Scalar             The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer         The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename _Scalar, typename _Vectorizer>
class USQOneElectronOperator:
    public SpinResolvedBase<USQOneElectronOperatorComponent<_Scalar, _Vectorizer>, USQOneElectronOperator<_Scalar, _Vectorizer>>,
    public SpinResolvedBasisTransformable<USQOneElectronOperator<_Scalar, _Vectorizer>>,
    public VectorSpaceArithmetic<USQOneElectronOperator<_Scalar, _Vectorizer>, _Scalar> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of 'this'.
    using Self = USQOneElectronOperator<Scalar, Vectorizer>;

    // The spinor tag corresponding to a `USQOneElectronOperator`.
    using SpinorTag = UnrestrictedSpinOrbitalTag;

    // The type of transformation matrix that is naturally related to a `USQOneElectronOperator`.
    using TM = UTransformationMatrix<Scalar>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<USQOneElectronOperatorComponent<Scalar, Vectorizer>, USQOneElectronOperator<Scalar, Vectorizer>>::SpinResolvedBase;


    /**
     *  Create an `USQOneElectronOperator` from all the matrix representations of its components.
     * 
     *  @param fs_a         All the matrix representations of the components of the alpha-part of the unrestricted one-electron operator.
     *  @param fs_b         All the matrix representations of the components of the beta-part of the unrestricted one-electron operator.
     */
    USQOneElectronOperator(const std::array<QCMatrix<Scalar>, Components>& fs_a, const std::array<QCMatrix<Scalar>, Components>& fs_b) :

        // Encapsulate the array of matrix representations in the alpha- and beta- operator components, and put them together to form the `USQOneElectronOperator`.
        SpinResolvedBase<USQOneElectronOperatorComponent<Scalar, Vectorizer>, USQOneElectronOperator<Scalar, Vectorizer>>(
            USQOneElectronOperatorComponent<Scalar, Vectorizer> {fs_a}, USQOneElectronOperatorComponent<Scalar, Vectorizer> {fs_b}) {

        // Check if the given matrix representations have the same dimensions.
        const auto dimension_of_first_a = fs_a[0].numberOfOrbitals();
        const auto dimension_of_first_b = fs_b[0].numberOfOrbitals();

        for (size_t i = 1; i < Components; i++) {
            const auto dimension_of_ith_a = fs_a[i].numberOfOrbitals();
            const auto dimension_of_ith_b = fs_b[i].numberOfOrbitals();

            if ((dimension_of_first_a != dimension_of_ith_a) || (dimension_of_first_b != dimension_of_ith_b)) {
                throw std::invalid_argument("USQOneElectronOperator(const std::array<QCMatrix<Scalar>, Components>&, const std::array<QCMatrix<Scalar>, components>&): The given matrix representations do not have the same dimensions for either the alpha or beta component.");
            }
        }
    }


    /**
     *  A constructor for ScalarUSQOneElectronOperators that doesn't require the argument to be a vector of just one element.
     *
     *  @param f_a          The matrix representation of the alpha-part of the unrestricted one-electron operator.
     *  @param f_b          The matrix representation of the beta-part of the unrestricted one-electron operator.
     *
     *  @note This constructor is only available for ScalarUSQOneElectronOperators (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415).
     */
    template <typename Z = Vectorizer>
    USQOneElectronOperator(const QCMatrix<Scalar>& f_a, const QCMatrix<Scalar>& f_b,
                           typename std::enable_if<std::is_same<Z, ScalarVectorizer>::value>::type* = 0) :
        USQOneElectronOperator(std::array<QCMatrix<Scalar>, 1> {f_a}, std::array<QCMatrix<Scalar>, 1> {f_b}) {}


    /**
     *  The default constructor.
     */
    USQOneElectronOperator() :
        USQOneElectronOperator(USQOneElectronOperator<Scalar, Vectorizer>::Zero(0))  // The dimensions of the representations are zero.
    {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Construct an `USQOneElectronOperator` with parameters that are zero, for both the alpha and beta spin components.
     *
     *  @param dim          The dimension of the matrix representation of the parameters, equal for the alpha and beta component.
     */
    static USQOneElectronOperator<Scalar, Vectorizer> Zero(const size_t dim) {

        const auto zero_component = USQOneElectronOperatorComponent<Scalar, Vectorizer>::Zero(dim);
        return USQOneElectronOperator<Scalar, Vectorizer> {zero_component, zero_component};
    }


    /*
     *  MARK: Parameter access
     */

    /**
     *  Access a component of this operator.
     * 
     *  @param indices      A set of coordinates that accesses this operator.
     * 
     *  @return The component of this operator that corresponds to the given coordinate indices.
     */
    template <typename... Indices>
    USQOneElectronOperator<Scalar, ScalarVectorizer> operator()(const Indices&... indices) const {

        return USQOneElectronOperator<Scalar, ScalarVectorizer> {this->alpha().parameters(indices), this->beta().parameters(indices)};
    }


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
    StorageArray<Scalar, Vectorizer> calculateExpectationValue(const SpinResolved1DM<Scalar>& D) const {

        // Calculate the sum of the alpha- and beta-contributions.
        return this->alpha().calculateExpectationValue(D.alpha()) + this->beta().calculateExpectationValue(D.beta());
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of alpha- and beta-orbitals.
     */
    size_t numberOfOrbitals() const {

        return this->alpha().numberOfOrbitals() + this->beta().numberOfOrbitals();
    }


    /**
     *  @param sigma                Alpha or beta.
     *
     *  @return The number of orbitals for the given spin component.
     */
    size_t numberOfOrbitals(const Spin sigma) const { return this->component(sigma).numberOfOrbitals(); }


    /**
     *  Apply the basis transformation and return the result.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed object.
     */
    // Self transformed(const TM& transformation_matrix) const {

    //     auto result = *this;
    // }


    /*
     *  MARK: Conforming to VectorSpaceArithmetic
     */

    /**
     *  Addition-assignment, to be implemented in derived classes.
     */
    Self& operator+=(const Self& rhs) {

        // Add the alpha-components and the beta-components.
        this->alpha() += rhs.alpha();
        this->beta() += rhs.beta();

        return *this;
    }


    /**
     *  Scalar multiplication-assignment, to be implemented in derived classes.
     */
    Self& operator*=(const Scalar& a) {

        // Multiply the alpha- and beta-components with the scalar.
        this->alpha() *= a;
        this->beta() *= a;

        return *this;
    }
};


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
struct OperatorTraits<USQSQOneElectronOperator<_Scalar, _Vectorizer, _DerivedOperator>> {

    // The scalar type used for a single parameter/matrix element/integral: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator that derives from `SimpleSQOneElectronOperator`, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = _DerivedOperator;
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like USQOneElectronOperator, i.e. with scalar-like access.
template <typename Scalar>
using ScalarRSQOneElectronOperator = USQOneElectronOperator<Scalar, ScalarVectorizer>;

// A vector-like USQOneElectronOperator, i.e. with vector-like access.
template <typename Scalar>
using VectorRSQOneElectronOperator = USQOneElectronOperator<Scalar, VectorVectorizer>;

// A matrix-like USQOneElectronOperator, i.e. with matrix-like access.
template <typename Scalar>
using MatrixRSQOneElectronOperator = USQOneElectronOperator<Scalar, MatrixVectorizer>;

// A tensor-like USQOneElectronOperator, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorRSQOneElectronOperator = USQOneElectronOperator<Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information (traits) on `USQOneElectronOperator` that is otherwise not accessible through a public class alias.
 * 
 *  @tparam Scalar          The scalar type used for a single parameter: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename Scalar, typename Vectorizer>
struct OperatorTraits<USQOneElectronOperator<Scalar, Vectorizer>> {

    // A type that corresponds to the scalar version of the associated unrestricted one-electron operator type.
    using ScalarOperator = ScalarUSQOneElectronOperator<Scalar>;

    // The type of transformation matrix that is naturally associated to an unrestricted one-electron operator.
    using TM = UTransformationMatrix<Scalar>;

    // The type of the one-particle density matrix that is naturally associated an unrestricted one-electron operator.
    using OneDM = SpinResolved1DM<Scalar>;

    // The type of the two-particle density matrix that is naturally associated an unrestricted one-electron operator.
    using TwoDM = SpinResolved2DM<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar, typename Vectorizer>
struct BasisTransformableTraits<USQOneElectronOperator<Scalar, Vectorizer>> {

    // The type of transformation matrix that is naturally related to a `USQOneElectronOperator`.
    using TM = UTransformationMatrix<Scalar>;
};


}  // namespace GQCP
