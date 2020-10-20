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


#include "Basis/Transformations/UTransformationMatrix.hpp"
#include "DensityMatrix/G1DM.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "DensityMatrix/SpinDensity1DM.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"
#include "Mathematical/Functions/VectorSpaceArithmetic.hpp"
#include "QuantumChemical/SpinResolved.hpp"


namespace GQCP {


/**
 *  A type that encapsulates alpha-alpha and beta-beta (spin-resolved) density matrices.
 *
 *  @tparam _Scalar             The scalar type of one of the density matrix elements: real or complex.
 */
template <typename _Scalar>
class SpinResolved1DM:
    public SpinResolvedBase<SpinResolved1DMComponent<_Scalar>, SpinResolved1DM<_Scalar>>,
    public BasisTransformable<SpinResolved1DM<_Scalar>>,
    public VectorSpaceArithmetic<SpinResolved1DM<_Scalar>, _Scalar> {
public:
    // The scalar type of one of the density matrix elements: real or complex.
    using Scalar = _Scalar;

    // The type of the transformation matrix that is naturally related to SpinResolved1DM.
    using TM = UTransformationMatrix<Scalar>;

    // The type of 'this'.
    using Self = SpinResolved1DM<Scalar>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<SpinResolved1DMComponent<_Scalar>, SpinResolved1DM<_Scalar>>::SpinResolvedBase;


    /*
     *  MARK: Constructors
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
    size_t numberOfOrbitals(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->alpha().numberOfOrbitals();
        }
        case Spin::beta: {
            return this->beta().numberOfOrbitals();
        }
        }
    }


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
    Orbital1DM<Scalar> spinSummed() const {
        return this->alpha() + this->beta();
    }


    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the result.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed object.
     */
    SpinResolved1DM<Scalar> transformed(const UTransformationMatrix<Scalar>& transformation_matrix) const override {

        auto result = *this;

        // Transform the components with the components of the transformation matrix.
        result.alpha().transform(transformation_matrix.alpha());
        result.beta().transform(transformation_matrix.beta());

        return result;
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
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<SpinResolved1DM<Scalar>> {

    // The type of the transformation matrix that is naturally related to SpinResolved1DM.
    using TM = UTransformationMatrix<Scalar>;
};


}  // namespace GQCP
