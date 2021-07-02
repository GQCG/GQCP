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


#include "Basis/Transformations/GTransformation.hpp"
#include "Basis/Transformations/JacobiRotation.hpp"
#include "DensityMatrix/G1DM.hpp"
#include "DensityMatrix/G2DM.hpp"
#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/PureUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/SimpleSQTwoElectronOperator.hpp"
#include "QuantumChemical/spinor_tags.hpp"


namespace GQCP {


/**
 *  A general(ized) two-electron operator, which is suited for expressing spin-dependent two-electron operators.
 *
 *  @tparam _Scalar                 The scalar type used for a single parameter/matrix element: real or complex.
 *  @tparam _Vectorizer             The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
 */
template <typename _Scalar, typename _Vectorizer>
class GSQTwoElectronOperator:
    public SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, GSQTwoElectronOperator<_Scalar, _Vectorizer>> {
public:
    // The scalar type used for a single parameter/matrix element: real or complex.
    using Scalar = _Scalar;

    //The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
    using Vectorizer = _Vectorizer;

    // The spinor tag corresponding to a `GSQTwoElectronOperator`.
    using SpinorTag = GeneralSpinorTag;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SimpleSQOneElectronOperator`'s constructors.
    using SimpleSQTwoElectronOperator<_Scalar, _Vectorizer, GSQTwoElectronOperator<_Scalar, _Vectorizer>>::SimpleSQTwoElectronOperator;


    /*
     *  MARK: Named constructors
     */

    /**
     *  Construct a `GSQTwoElectronOperator` from a `PureUSQTwoElectronOperatorComponent`.
     *
     *  @param g_component          The pure component of an unrestricted two-electron operator that should be converted.
     */
    static GSQTwoElectronOperator<Scalar, Vectorizer> FromUnrestrictedComponent(const PureUSQTwoElectronOperatorComponent<Scalar, Vectorizer>& g_component) {

        // We can just wrap the one-electron integrals into the correct class.
        const StorageArray<SquareRankFourTensor<Scalar>, Vectorizer> array {g_component.allParameters(), g_component.vectorizer()};
        return GSQTwoElectronOperator<Scalar, Vectorizer> {array};
    }
};


/*
 *  MARK: Convenience aliases
 */

// A scalar-like GSQTwoElectronOperator, i.e. with scalar-like access.
template <typename Scalar>
using ScalarGSQTwoElectronOperator = GSQTwoElectronOperator<Scalar, ScalarVectorizer>;

// A vector-like GSQTwoElectronOperator, i.e. with vector-like access.
template <typename Scalar>
using VectorGSQTwoElectronOperator = GSQTwoElectronOperator<Scalar, VectorVectorizer>;

// A matrix-like GSQTwoElectronOperator, i.e. with matrix-like access.
template <typename Scalar>
using MatrixGSQTwoElectronOperator = GSQTwoElectronOperator<Scalar, MatrixVectorizer>;

// A tensor-like GSQTwoElectronOperator, i.e. with tensor-like access.
template <typename Scalar, size_t N>
using TensorGSQTwoElectronOperator = GSQTwoElectronOperator<Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information (traits) on `GSQTwoElectronOperator` that is otherwise not accessible through a public class alias.
 *
 *  @tparam Scalar          The scalar type used for a single parameter/matrix element: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
 */
template <typename Scalar, typename Vectorizer>
struct OperatorTraits<GSQTwoElectronOperator<Scalar, Vectorizer>> {

    // A type that corresponds to the scalar version of the associated general(ized) two-electron operator type.
    using ScalarOperator = ScalarGSQTwoElectronOperator<Scalar>;

    // The type of one-electron operator that is naturally related to a `GSQTwoElectronOperator`.
    using SQOneElectronOperator = GSQOneElectronOperator<Scalar, Vectorizer>;

    // The type of transformation that is naturally associated to a `GSQTwoElectronOperator`.
    using Transformation = GTransformation<Scalar>;

    // The type of density matrix that is naturally associated to a `GSQTwoElectronOperator`.
    using OneDM = G1DM<Scalar>;

    // The type of density matrix that is naturally associated to a `GSQTwoElectronOperator`.
    using TwoDM = G2DM<Scalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 *
 *  @tparam Scalar          The scalar type used for a single parameter/matrix element: real or complex.
 *  @tparam Vectorizer      The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
 */
template <typename Scalar, typename Vectorizer>
struct BasisTransformableTraits<GSQTwoElectronOperator<Scalar, Vectorizer>> {

    // The type of transformation matrix that is naturally associated to a `GSQTwoElectronOperator`.
    using Transformation = GTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar, typename Vectorizer>
struct JacobiRotatableTraits<GSQTwoElectronOperator<Scalar, Vectorizer>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


/*
 *  MARK: One-electron operator products
 */

/**
 *  A type that encapsulates the matrix elements of the product of two scalar, generalized one-electron operators.
 *
 *  @tparam _Scalar         The scalar type of one of the matrix elements: real or complex.
 */
template <typename _Scalar>
class ScalarGSQOneElectronOperatorProduct:
    public VectorSpaceArithmetic<ScalarGSQOneElectronOperatorProduct<_Scalar>, _Scalar>,
    public BasisTransformable<ScalarGSQOneElectronOperatorProduct<_Scalar>>,
    public JacobiRotatable<ScalarGSQOneElectronOperatorProduct<_Scalar>> {
public:
    // The scalar type of one of the matrix elements: real or complex.
    using Scalar = _Scalar;

    // The type of 'this'.
    using Self = ScalarGSQOneElectronOperatorProduct<Scalar>;


private:
    // The one-electron part of the product.
    ScalarGSQOneElectronOperator<Scalar> o;

    // The two-electron part of the product. Since it is expressed as a `SQTwoElectronOperator`, it implicitly assumes the incorporation of a factor 2.
    ScalarGSQTwoElectronOperator<Scalar> t;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create the representation of the product of two scalar, generalized one-electron operators by its one- and two-electron constituents.
     *
     *  @param o            The one-electron part of the product.
     *  @param t            The two-electron part of the product. Since it is expressed as a `SQTwoElectronOperator`, it implicitly assumes the incorporation of a factor 2.
     */
    ScalarGSQOneElectronOperatorProduct(const ScalarGSQOneElectronOperator<Scalar>& o, const ScalarGSQTwoElectronOperator<Scalar>& t) :
        o {o},
        t {t} {}


    /*
     *  MARK: Access
     */

    /**
     *  @return A read-only reference to the one-electron part of the product.
     */
    const ScalarGSQOneElectronOperator<Scalar>& oneElectron() const { return this->o; }

    /**
     *  @return A read-only reference to the two-electron part of the product.
     */
    const ScalarGSQTwoElectronOperator<Scalar>& twoElectron() const { return this->t; }

    /**
     *  @return A writable reference to the one-electron part of the product.
     */
    ScalarGSQOneElectronOperator<Scalar>& oneElectron() { return this->o; }

    /**
     *  @return A writable reference to the two-electron part of the product.
     */
    ScalarGSQTwoElectronOperator<Scalar>& twoElectron() { return this->t; }


    /*
     *  MARK: Conforming to `VectorSpaceArithmetic`
     */

    /**
     *  Addition-assignment.
     */
    ScalarGSQOneElectronOperatorProduct& operator+=(const ScalarGSQOneElectronOperatorProduct& rhs) override {

        // Add the one- and two-electron operator parts.
        this->oneElectron() += rhs.oneElectron();
        this->twoElectron() += rhs.twoElectron();

        return *this;
    }


    /**
     *  Scalar multiplication-assignment.
     */
    ScalarGSQOneElectronOperatorProduct& operator*=(const Scalar& a) override {

        // Multiply the one- and two-electron operator parts with the scalar.
        this->oneElectron() *= a;
        this->twoElectron() *= a;

        return *this;
    }


    /*
     *  MARK: Conforming to BasisTransformable
     */

    /**
     *  Apply the basis transformation and return the resulting Hamiltonian.
     *
     *  @param T            The basis transformation.
     *
     *  @return The basis-transformed Hamiltonian.
     */
    Self transformed(const GTransformation<Scalar>& T) const override {

        auto result = *this;

        // Transform the one- and two-electron contributions.
        result.oneElectron().transform(T);
        result.twoElectron().transform(T);

        return result;
    }


    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<Self>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<Self>::rotated;


    /*
     *  MARK: Conforming to JacobiRotatable
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     *
     *  @param jacobi_rotation          The Jacobi rotation.
     *
     *  @return The Jacobi-rotated object.
     */
    Self rotated(const JacobiRotation& jacobi_rotation) const override {

        auto result = *this;

        // Transform the total one- and two-electron contributions.
        result.oneElectron().rotate(jacobi_rotation);
        result.twoElectron().rotate(jacobi_rotation);

        return result;
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<Self>::rotate;


    /*
     *  MARK: Expectation value
     */

    /**
     *  Calculate the expectation value of this one-electron operator product.
     *
     *  @param D            The 1-DM.
     *  @param d            The 2-DM.
     *
     *  @return The expectation value of this one-electron operator product.
     */
    Scalar calculateExpectationValue(const G1DM<Scalar>& D, const G2DM<Scalar>& d) const {

        // We can access the expectation value of a scalar operator using an empty call.
        return this->oneElectron().calculateExpectationValue(D)() + this->twoElectron().calculateExpectationValue(d)();
    }
};


/*
 *  MARK: `BasisTransformableTraits` for `ScalarGSQOneElectronOperatorProduct`.
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 *
 *  @tparam Scalar          The scalar type of one of the matrix elements: real or complex.
 */
template <typename Scalar>
struct BasisTransformableTraits<ScalarGSQOneElectronOperatorProduct<Scalar>> {

    // The type of transformation matrix that is naturally associated to a `ScalarGSQOneElectronOperatorProduct`.
    using Transformation = GTransformation<Scalar>;
};


/*
 *  MARK: `JacobiRotatableTraits` for `ScalarGSQOneElectronOperatorProduct`.
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 *
 *  @tparam Scalar           The scalar type of one of the matrix elements: real or complex.
 */
template <typename Scalar>
struct JacobiRotatableTraits<ScalarGSQOneElectronOperatorProduct<Scalar>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


/**
 *  Multiply one scalar, generalized one-electron operator by another.
 *
 *  @tparam Scalar           The scalar type of one of the matrix elements: real or complex.
 *
 *  @note This function assumes that
 */
template <typename Scalar>
ScalarGSQOneElectronOperatorProduct<Scalar> operator*(const ScalarGSQOneElectronOperator<Scalar>& lhs, const ScalarGSQOneElectronOperator<Scalar>& rhs) {

    // Prepare some variables.
    const auto& L = lhs.parameters();
    const auto& R = rhs.parameters();

    // Determine the matrix elements of the one-electron part of the product.
    SquareMatrix<Scalar> O = L * R;
    ScalarGSQOneElectronOperator<Scalar> O_op {O};


    // Determine the matrix elements of the two-electron part of the product.
    const auto dim = lhs.parameters().dimension();
    SquareRankFourTensor<Scalar> T = SquareRankFourTensor<Scalar>::Zero(dim);
    for (size_t p = 0; p < dim; p++) {
        for (size_t q = 0; q < dim; q++) {
            for (size_t r = 0; r < dim; r++) {
                for (size_t s = 0; s < dim; s++) {
                    T(p, q, r, s) = 2.0 * L(p, q) * R(r, s);  // Include the prefactor '2' because we're going to encapsulate these matrix elements with a `ScalarGSQTwoElectronOperator`, whose matrix elements should not embet the prefactor 0.5 for two-electron operators.
                }
            }
        }
    }
    ScalarGSQTwoElectronOperator<Scalar> T_op {T};


    return ScalarGSQOneElectronOperatorProduct<Scalar>(O_op, T_op);
}


}  // namespace GQCP
