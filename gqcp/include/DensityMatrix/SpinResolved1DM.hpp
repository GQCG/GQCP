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
#include "Basis/Transformations/SpinResolvedJacobiRotatable.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "DensityMatrix/G1DM.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "DensityMatrix/SpinDensity1DM.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"
#include "Mathematical/Functions/VectorSpaceArithmetic.hpp"
#include "QuantumChemical/SpinResolved.hpp"


namespace GQCP {


/**
 *  A type that encapsulates alpha and beta (spin-resolved) density matrices.
 *
 *  @tparam _Scalar             The scalar type of one of the density matrix elements: real or complex.
 */
template <typename _Scalar>
class SpinResolved1DM:
    public SpinResolvedBase<SpinResolved1DMComponent<_Scalar>, SpinResolved1DM<_Scalar>>,
    public SpinResolvedBasisTransformable<SpinResolved1DM<_Scalar>>,
    public SpinResolvedJacobiRotatable<SpinResolved1DM<_Scalar>>,
    public VectorSpaceArithmetic<SpinResolved1DM<_Scalar>, _Scalar> {
public:
    // The scalar type of one of the density matrix elements: real or complex.
    using Scalar = _Scalar;

    // The type of the transformation that is naturally related to a `SpinResolved1DM`.
    using Transformation = UTransformation<Scalar>;

    // The type of 'this'.
    using Self = SpinResolved1DM<Scalar>;

    // The type component this spin resolved object is made of.
    using ComponentType = typename SpinResolvedBase<SpinResolved1DMComponent<Scalar>, Self>::Of;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<SpinResolved1DMComponent<_Scalar>, SpinResolved1DM<_Scalar>>::SpinResolvedBase;


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create a spin-resolved 1-DM from an `Orbital1DM`, attributing half of the orbital 1-DM to each of the spin components.
     * 
     *  @param D            The orbital 1-DM.
     * 
     *  @return A spin-resolved 1-DM.
     */
    static SpinResolved1DM<Scalar> FromOrbital1DM(const Orbital1DM<Scalar>& D) {

        const SpinResolved1DMComponent<Scalar> D_half = D / 2;
        return SpinResolved1DM<Scalar>(D_half, D_half);
    }


    /*
     *  MARK: Conversions
     */

    /**
     *  @return This spin-resolved 1-DM as a generalized 1-DM (`G1DM`).
     * 
     *  @note We cannot implement this as a named constructor on `G1DM` because we require `norm` to be implemented on `SpinResolved1DM` and that internally uses a `G1DM`, hence we have to avoid the circular dependency.
     */
    G1DM<Scalar> generalized() const {

        // Determine the dimensions of the generalized, spin-blocked 1-DM.
        const auto K_alpha = this->alpha().numberOfOrbitals();
        const auto K_beta = this->beta().numberOfOrbitals();
        const auto M = K_alpha + K_beta;

        // The generalized 1-DM contains the alpha part in the top-left corner, and the beta part in the bottom right corner.
        G1DM<Scalar> D_generalized = G1DM<Scalar>::Zero(M);
        D_generalized.topLeftCorner(K_alpha, K_alpha) = this->alpha();
        D_generalized.bottomRightCorner(K_beta, K_beta) = this->beta();

        return D_generalized;
    }


    /*
     *  MARK: General information
     */

    /**
     * @return The Frobenius norm of this spin-resolved 1-DM.
     */
    double norm() const { return this->generalized().norm(); }

    /**
     *  @param sigma            Alpha or beta.
     * 
     *  @return The number of orbitals (spinors or spin-orbitals, depending on the context) that correspond to the given spin.
     */
    size_t numberOfOrbitals(const Spin sigma) const { return this->component(sigma).numberOfOrbitals(); }


    /*
     *  MARK: Spin-related operations
     */

    /**
     *  @return The spin-density matrix, i.e. the difference between the alpha and beta 1-DM.
     */
    SpinDensity1DM<Scalar> spinDensity() const {
        return this->alpha() - this->beta();
    }

    /**
     *  @return The orbital density matrix, i.e. the sum of the alpha and beta 1-DM.
     */
    Orbital1DM<Scalar> orbitalDensity() const {
        return this->alpha() + this->beta();
    }


    /*
     *  MARK: Conforming to `VectorSpaceArithmetic`
     */

    /**
     *  Addition-assignment.
     */
    Self& operator+=(const Self& rhs) override {

        // For addition, the alpha- and beta-parts should be added component-wise.
        this->alpha() += rhs.alpha();
        this->beta() += rhs.beta();

        return *this;
    }

    /**
     *  Scalar multiplication-assignment
     */
    Self& operator*=(const Scalar& a) override {

        // For scalar multiplication, the alpha- and beta-parts should be multiplied with the scalar.
        this->alpha() *= a;
        this->beta() *= a;

        return *this;
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
struct BasisTransformableTraits<SpinResolved1DM<Scalar>> {

    // The type of the transformation that is naturally related to a `SpinResolved1DM`.
    using Transformation = UTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<SpinResolved1DM<Scalar>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = UJacobiRotation;
};


}  // namespace GQCP
