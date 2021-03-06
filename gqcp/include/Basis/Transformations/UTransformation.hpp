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
#include "Basis/Transformations/SpinResolvedBasisTransformable.hpp"
#include "Basis/Transformations/SpinResolvedJacobiRotatable.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/*
 *  MARK: UTransformation implementation
 */

/**
 *  A type that encapsulates transformation matrices for the alpha- and beta-parts of spin-orbital bases.
 * 
 *  @tparam _Scalar         The scalar type used for a transformation coefficient: real or complex.
 */
template <typename _Scalar>
class UTransformation:
    public SpinResolvedBase<UTransformationComponent<_Scalar>, UTransformation<_Scalar>>,
    public SpinResolvedBasisTransformable<UTransformation<_Scalar>>,
    public SpinResolvedJacobiRotatable<UTransformation<_Scalar>> {
public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;

    // The type of the transformation for which the basis transformation should be defined. A transformation should naturally be transformable with itself.
    using Transformation = UTransformation<Scalar>;

    // The type of 'this'.
    using Self = UTransformation<Scalar>;

    // The type component this spin resolved object is made of.
    using ComponentType = typename SpinResolvedBase<UTransformationComponent<Scalar>, Self>::Of;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<UTransformationComponent<Scalar>, UTransformation<Scalar>>::SpinResolvedBase;


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create an UTransformation from an RTransformation, leading to transformations for both spin components that are equal.
     * 
     *  @param T                The transformation that is equal for both the alpha and the beta spin-orbitals.
     * 
     *  @return A `UTransformation` corresponding to the given `RTransformation`.
     */
    static UTransformation<Scalar> FromRestricted(const RTransformation<Scalar>& T) {

        // Wrap the restricted transformation in the correct class and return the result.
        const UTransformationComponent<Scalar> T_component {T.matrix()};
        return UTransformation<Scalar>::FromEqual(T_component);
    }


    /**
     *  Create an identity UTransformation.
     * 
     *  @param dim_alpha            The number of alpha spin-orbitals.
     *  @param dim_beta             The number of beta spin-orbitals.
     * 
     *  @return An identity UTransformation.
     */
    static UTransformation<Scalar> Identity(const size_t dim_alpha, const size_t dim_beta) {

        const auto T_alpha = UTransformationComponent<Scalar>::Identity(dim_alpha);
        const auto T_beta = UTransformationComponent<Scalar>::Identity(dim_beta);

        return UTransformation<Scalar> {T_alpha, T_beta};
    }


    /**
     *  Create an identity UTransformation.
     * 
     *  @param dim              The dimension of the alpha and beta spin-orbitals.
     * 
     *  @return An identity UTransformation.
     */
    static UTransformation<Scalar> Identity(const size_t dim) { return UTransformation<Scalar>::Identity(dim, dim); }


    /**
     *  Create a random `UTransformation`.
     * 
     *  @param dim              The number of alpha or beta spin-orbitals (equal).
     * 
     *  @return A random `UTransformation`.
     */
    static UTransformation<Scalar> Random(const size_t dim) {

        const auto T_alpha = UTransformationComponent<Scalar>::Random(dim);
        const auto T_beta = UTransformationComponent<Scalar>::Random(dim);

        return UTransformation<Scalar> {T_alpha, T_beta};
    }


    /**
     *  Create a random unitary `UTransformation`.
     * 
     *  @param dim              The number of alpha or beta spin-orbitals (equal).
     * 
     *  @return A random unitary `UTransformation`.
     */
    static UTransformation<Scalar> RandomUnitary(const size_t dim) {

        const auto T_alpha = UTransformationComponent<Scalar>::RandomUnitary(dim);
        const auto T_beta = UTransformationComponent<Scalar>::RandomUnitary(dim);

        return UTransformation<Scalar> {T_alpha, T_beta};
    }


    /*
     *  MARK: Linear algebraic operations
     */

    /**
     *  @param threshold        The threshold used for checking unitarity.
     * 
     *  @return If this transformation matrix is considered to be unitary, within the given treshold.
     */
    bool isUnitary(const double threshold) const { return this->alpha().isUnitary(threshold) && this->beta().isUnitary(threshold); }


    /**
     *  @return The inverse of this transformation.
     */
    UTransformation<Scalar> inverse() const {

        auto result = *this;

        // The inverse of this compound transformation is the transformation where the alpha- and beta-transformations are inverted.
        result.alpha() = this->alpha().inverse();
        result.beta() = this->beta().inverse();

        return result;
    }


    /**
     *  MARK: Enabling basis transformations
     */

    // Since `rotate` and `rotated` are both defined in `SpinResolvedBasisTransformable` and `SpinResolvedJacobiRotatable`, we have to explicitly enable these methods here.

    // Allow the `rotate` method from `SpinResolvedBasisTransformable`, since there's also a `rotate` from `SpinResolvedJacobiRotatable`.
    using SpinResolvedBasisTransformable<Self>::rotate;

    // Allow the `rotated` method from `SpinResolvedBasisTransformable`, since there's also a `rotated` from `SpinResolvedJacobiRotatable`.
    using SpinResolvedBasisTransformable<Self>::rotated;

    // Allow the `rotate` method from `SpinResolvedJacobiRotatable`, since there's also a `rotate` from `SpinResolvedBasisTransformable`.
    using SpinResolvedJacobiRotatable<Self>::rotate;

    // Allow the `rotated` method from `SpinResolvedJacobiRotatable`, since there's also a `rotated` from `SpinResolvedBasisTransformable`.
    using SpinResolvedJacobiRotatable<Self>::rotated;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<UTransformation<Scalar>> {

    // The type of the transformation for which the basis transformation should be defined. A transformation matrix should naturally be transformable with itself.
    using Transformation = UTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<UTransformation<Scalar>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = UJacobiRotation;
};


}  // namespace GQCP
