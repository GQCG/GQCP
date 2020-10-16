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


#include "DensityMatrix/G1DM.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "DensityMatrix/SpinDensity1DM.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"
#include "QuantumChemical/SpinResolved.hpp"


namespace GQCP {


/**
 *  A type that encapsulates alpha-alpha and beta-beta (spin-resolved) density matrices.
 *
 *  @tparam _Scalar             The scalar type of one of the density matrix elements: real or complex.
 */
template <typename _Scalar>
class SpinResolved1DM:
    public SpinResolvedBase<SpinResolved1DMComponent<_Scalar>, SpinResolved1DM<_Scalar>> {
public:
    // The scalar type of one of the density matrix elements: real or complex.
    using Scalar = _Scalar;


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
     *  @param T_a          transformation matrix for the alpha component of the spin resolved 1-DM
     *  @param T_b          transformation matrix for the beta component of the spin resolved 1-DM
     * 
     *  @return the transformed spin resolved density matrix, with each component transformed seperately t a different basis.
     */
    SpinResolved1DM<Scalar> transformed(const TransformationMatrix<double>& T_a, const TransformationMatrix<double>& T_b) const {
        OneDM<Scalar> D_a_transformed = T_a.conjugate() * this->alpha() * T_a.transpose();
        OneDM<Scalar> D_b_transformed = T_b.conjugate() * this->beta() * T_b.transpose();

        return SpinResolved1DM<Scalar> {D_a_transformed, D_b_transformed};
    }

    /**
     *  @param T          transformation matrix for the alpha and beta component of the spin resolved 1-DM
     * 
     *  @return the transformed spin resolved density matrix, with each component transformed to the same basis.
     */
    SpinResolved1DM<Scalar> transform(const TransformationMatrix<double>& T) const {
        return this->transform(T, T);
    }
};


/*
*  OPERATORS
*/

/**
*  Add two spin resolved density matrices by adding their parameters. The two alphas are added together and the two betas are added together. 
* 
*  @tparam LHSScalar           the scalar type of the left-hand side
*  @tparam RHSScalar           the scalar type of the right-hand side
* 
*  @param lhs                  the left-hand side
*  @param rhs                  the right-hand side
*/
template <typename LHSScalar, typename RHSScalar>
auto operator+(const SpinResolved1DM<LHSScalar>& lhs, const SpinResolved1DM<RHSScalar>& rhs) -> SpinResolved1DM<sum_t<LHSScalar, RHSScalar>> {

    using ResultScalar = sum_t<LHSScalar, RHSScalar>;

    auto D_sum_a = lhs.alpha();
    auto D_sum_b = lhs.beta();

    D_sum_a += rhs.alpha();
    D_sum_b += rhs.beta();

    return SpinResolved1DM<ResultScalar>(D_sum_a, D_sum_b);
}

/**
 *  Multiply a one-electron operator with a scalar. The alpha and beta components are multiplied with the same scalar.
 * 
 *  @tparam Scalar                            the scalar type of the scalar
 *  @tparam DMScalar                          the scalar type of the spin resolved one DM
 * 
 *  @tparam scalar                            the scalar of the scalar multiplication
 *  @tparam spin_resolved_DM                  the spin resolved one DM
 */
template <typename Scalar, typename DMScalar>
auto operator*(const Scalar& scalar, const SpinResolved1DM<DMScalar>& spin_resolved_DM) -> SpinResolved1DM<product_t<Scalar, DMScalar>> {

    using ResultScalar = product_t<Scalar, DMScalar>;

    auto D_a = spin_resolved_DM.alpha();
    auto D_b = spin_resolved_DM.beta();

    D_a *= scalar;
    D_b *= scalar;

    return SpinResolved1DM<ResultScalar>(D_a, D_b);
}


/**
 *  Negate a one-electron operator
 * 
 *  @tparam Scalar              the scalar type of the spin resolved one DM
 * 
 *  @param spin_resolved_DM                   the spin resolved density matrix
 */
template <typename Scalar>
SpinResolved1DM<Scalar> operator-(const SpinResolved1DM<Scalar>& spin_resolved_DM) {

    return (-1.0) * spin_resolved_DM;  // negation is scalar multiplication with (-1.0)
}


/**
    *  Subtract two spin resolved density matricessubtracting their parameters
    * 
    *  @tparam LHSScalar           the scalar type of the left-hand side
    *  @tparam RHSScalar           the scalar type of the right-hand side
    * 
    *  @param lhs                  the left-hand side
    *  @param rhs                  the right-hand side
*/
template <typename LHSScalar, typename RHSScalar>
auto operator-(const SpinResolved1DM<LHSScalar>& lhs, const SpinResolved1DM<RHSScalar>& rhs) -> SpinResolved1DM<sum_t<LHSScalar, RHSScalar>> {

    return lhs + (-rhs);
}


}  // namespace GQCP
