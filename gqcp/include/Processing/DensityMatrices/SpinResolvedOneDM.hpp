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


#include "Basis/SpinorBasis/Spin.hpp"
#include "Processing/DensityMatrices/OneDM.hpp"


namespace GQCP {


/**
 *  A type that encapsulates alpha-alpha and beta-beta spin-resolved density matrices.
 *
 *  @tparam _Scalar             the scalar type of one of the elements
 */
template <typename _Scalar>
class SpinResolvedOneDM {
public:
    using Scalar = _Scalar;


private:
    OneDM<Scalar> D_aa;  // the alpha-alpha 1-DM
    OneDM<Scalar> D_bb;  // the beta-beta 1-DM


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Create a SpinResolvedOneDM from its members.
     *
     *  @param D_aa             the alpha-alpha 1-DM
     *  @param D_bb             the beta-beta 1-DM
     */
    SpinResolvedOneDM(const OneDM<Scalar>& D_aa, const OneDM<Scalar>& D_bb) :
        D_aa {D_aa},
        D_bb {D_bb} {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Create a spin-resolved 1-DM as half of the total 1-DM.
     * 
     *  @param D            the spin-summed 1-DM
     * 
     *  @return a spin-resolved 1-DM
     */
    static SpinResolvedOneDM<Scalar> FromRestricted(const OneDM<Scalar>& D) {

        const auto D_half = D / 2;
        return SpinResolvedOneDM<Scalar>(D_half, D_half);
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the alpha-alpha part of the spin-resolved 1-DM
     */
    const OneDM<Scalar>& alpha() const { return this->D_aa; }

    /**
     *  @return the beta-beta part of the spin-resolved 1-DM
     */
    const OneDM<Scalar>& beta() const { return this->D_bb; }

    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the number of orbitals (spinors or spin-orbitals, depending on the context) that correspond to the given spin
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

    /**
     *  @return the spin-density matrix, i.e. the difference between the alpha and beta 1-DM
     */
    OneDM<Scalar> spinDensity() const {
        return this->alpha() - this->beta();
    }

    /**
     *  @return the spin-summed density matrix, i.e. the sum of the alpha and beta 1-DM
     */
    OneDM<Scalar> spinSummed() const {
        return this->alpha() + this->beta();
    }

    /**
     *  @param T          transformation matrix for the alpha and beta component of the spin resolved 1-DM
     * 
     *  @return the transformed spin resolved density matrix, with each component transformed to the same basis.
     */
    SpinResolvedOneDM<Scalar> transform(const TransformationMatrix<double>& T) const {
        OneDM<Scalar> D_a_transformed = T.conjugate() * this->alpha() * T.transpose();
        OneDM<Scalar> D_b_transformed = T.conjugate() * this->beta() * T.transpose();

        return SpinResolvedOneDM<Scalar> {D_a_transformed, D_b_transformed};
    }

    /**
     *  @param T_a          transformation matrix for the alpha component of the spin resolved 1-DM
     *  @param T_b          transformation matrix for the beta component of the spin resolved 1-DM
     * 
     *  @return the transformed spin resolved density matrix, with each component transformed seperately.
     */
    SpinResolvedOneDM<Scalar> transform(const TransformationMatrix<double>& T_a, const TransformationMatrix<double>& T_b) const {
        OneDM<Scalar> D_a_transformed = T_a.conjugate() * this->alpha() * T_a.transpose();
        OneDM<Scalar> D_b_transformed = T_b.conjugate() * this->beta() * T_b.transpose();

        return SpinResolvedOneDM<Scalar> {D_a_transformed, D_b_transformed};
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
auto operator+(const SpinResolvedOneDM<LHSScalar>& lhs, const SpinResolvedOneDM<RHSScalar>& rhs) -> SpinResolvedOneDM<sum_t<LHSScalar, RHSScalar>> {

    using ResultScalar = sum_t<LHSScalar, RHSScalar>;

    auto D_sum_a = lhs.alpha();
    auto D_sum_b = lhs.beta();

    D_sum_a += rhs.alpha();
    D_sum_b += rhs.beta();

    return SpinResolvedOneDM<ResultScalar>(D_sum_a, D_sum_b);
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
auto operator*(const Scalar& scalar, const SpinResolvedOneDM<DMScalar>& spin_resolved_DM) -> SpinResolvedOneDM<product_t<Scalar, DMScalar>> {

    using ResultScalar = product_t<Scalar, DMScalar>;

    auto D_a = spin_resolved_DM.alpha();
    auto D_b = spin_resolved_DM.beta();

    D_a *= scalar;
    D_b *= scalar;

    return SpinResolvedOneDM<ResultScalar>(D_a, D_b);
}


/**
 *  Negate a one-electron operator
 * 
 *  @tparam Scalar              the scalar type of the spin resolved one DM
 * 
 *  @param spin_resolved_DM                   the spin resolved density matrix
 */
template <typename Scalar>
SpinResolvedOneDM<Scalar> operator-(const SpinResolvedOneDM<Scalar>& spin_resolved_DM) {

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
auto operator-(const SpinResolvedOneDM<LHSScalar>& lhs, const SpinResolvedOneDM<RHSScalar>& rhs) -> SpinResolvedOneDM<sum_t<LHSScalar, RHSScalar>> {

    return lhs + (-rhs);
}


}  // namespace GQCP
