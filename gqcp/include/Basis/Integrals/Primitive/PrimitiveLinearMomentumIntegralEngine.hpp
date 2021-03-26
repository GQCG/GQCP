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

#include "Basis/Integrals/Primitive/BaseVectorPrimitiveIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveOverlapIntegralEngine.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Operator/FirstQuantized/LinearMomentumOperator.hpp"
#include "Utilities/literals.hpp"


namespace GQCP {


/**
 *  A class that can calculate linear momentum integrals.
 * 
 *  @tparam _Shell              The type of shell that this integral engine is related to.
 */
template <typename _Shell>
class PrimitiveLinearMomentumIntegralEngine:
    public BaseVectorPrimitiveIntegralEngine {
public:
    // The type of shell that this integral engine is related to.
    using Shell = _Shell;

    // The type of primitive that underlies the type of shell.
    using Primitive = typename Shell::Primitive;

    // The scalar representation of a linear momentum integral.
    using IntegralScalar = product_t<LinearMomentumOperator::Scalar, typename Primitive::OutputType>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `BaseVectorPrimitiveIntegralEngine`'s constructors.
    using BaseVectorPrimitiveIntegralEngine::BaseVectorPrimitiveIntegralEngine;


    /*
     *  MARK: CartesianGTO integrals
     */

    /**
     *  Calculate the linear momentum integral (of the current component) over the two Cartesian GTOs.
     * 
     *  @param left             The left Cartesian GTO.
     *  @param right            The right Cartesian GTO.
     * 
     *  @return The linear momentum integral over the two given Cartesian GTOs.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, GTOShell>::value, IntegralScalar> calculate(const CartesianGTO& left, const CartesianGTO& right) {

        // Prepare some variables.
        const auto i = static_cast<int>(left.cartesianExponents().value(CartesianDirection::x));
        const auto k = static_cast<int>(left.cartesianExponents().value(CartesianDirection::y));
        const auto m = static_cast<int>(left.cartesianExponents().value(CartesianDirection::z));

        const auto j = static_cast<int>(right.cartesianExponents().value(CartesianDirection::x));
        const auto l = static_cast<int>(right.cartesianExponents().value(CartesianDirection::y));
        const auto n = static_cast<int>(right.cartesianExponents().value(CartesianDirection::z));

        const auto a = left.gaussianExponent();
        const auto b = right.gaussianExponent();

        const auto K_x = left.center()(CartesianDirection::x);
        const auto K_y = left.center()(CartesianDirection::y);
        const auto K_z = left.center()(CartesianDirection::z);

        const auto L_x = right.center()(CartesianDirection::x);
        const auto L_y = right.center()(CartesianDirection::y);
        const auto L_z = right.center()(CartesianDirection::z);

        PrimitiveOverlapIntegralEngine<GTOShell> S;


        // For the current component, the integral can be calculated as a product of three contributions.
        switch (this->component) {
        case CartesianDirection::x: {
            return this->calculate1D(a, K_x, i, b, L_x, j) * S.calculate1D(a, K_y, k, b, L_y, l) * S.calculate1D(a, K_z, m, b, L_z, n);
            break;
        }

        case CartesianDirection::y: {
            return S.calculate1D(a, K_x, i, b, L_x, j) * this->calculate1D(a, K_y, k, b, L_y, l) * S.calculate1D(a, K_z, m, b, L_z, n);
            break;
        }

        case CartesianDirection::z: {
            return S.calculate1D(a, K_x, i, b, L_x, j) * S.calculate1D(a, K_y, k, b, L_y, l) * this->calculate1D(a, K_z, m, b, L_z, n);
            break;
        }
        }
    }


    /**
     *  Calculate the linear momentum integral over two Cartesian GTO 1-D primitives.
     * 
     *  @param a                The Gaussian exponent of the left 1-D primitive.
     *  @param K                The (directional coordinate of the) center of the left 1-D primitive.
     *  @param i                The Cartesian exponent of the left 1-D primitive.
     *  @param b                The Gaussian exponent of the right 1-D primitive.
     *  @param L                The (directional coordinate of the) center of the right 1-D primitive.
     *  @param j                The Cartesian exponent of the right 1-D primitive.
     * 
     *  @return The linear momentum integral over the two given 1-D primitives.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, GTOShell>::value, IntegralScalar> calculate1D(const double a, const double K, const int i, const double b, const double L, const int j) {

        // The linear momentum integral is expressed entirely using overlap integrals.
        PrimitiveOverlapIntegralEngine<GTOShell> S;

        using namespace GQCP::literals;
        return 2.0 * 1.0_ii * b * S.calculate1D(a, K, i, b, L, j + 1) -
               1.0_ii * static_cast<double>(j) * S.calculate1D(a, K, i, b, L, j - 1);
    }


    /*
     *  MARK: London CartesianGTO integrals
     */

    /**
     *  Calculate the linear momentum integral over two London Cartesian GTOs.
     * 
     *  @param left             The left London Cartesian GTO.
     *  @param right            The right London Cartesian GTO.
     * 
     *  @return The linear momentum integral over the two given London Cartesian GTOs.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, IntegralScalar> calculate(const LondonCartesianGTO& left, const LondonCartesianGTO& right) {

        // Prepare some variables.
        const auto i = static_cast<int>(left.cartesianGTO().cartesianExponents().value(CartesianDirection::x));
        const auto k = static_cast<int>(left.cartesianGTO().cartesianExponents().value(CartesianDirection::y));
        const auto m = static_cast<int>(left.cartesianGTO().cartesianExponents().value(CartesianDirection::z));

        const auto j = static_cast<int>(right.cartesianGTO().cartesianExponents().value(CartesianDirection::x));
        const auto l = static_cast<int>(right.cartesianGTO().cartesianExponents().value(CartesianDirection::y));
        const auto n = static_cast<int>(right.cartesianGTO().cartesianExponents().value(CartesianDirection::z));

        const auto a = left.cartesianGTO().gaussianExponent();
        const auto b = right.cartesianGTO().gaussianExponent();

        const auto K_x = left.cartesianGTO().center()(CartesianDirection::x);
        const auto K_y = left.cartesianGTO().center()(CartesianDirection::y);
        const auto K_z = left.cartesianGTO().center()(CartesianDirection::z);

        const auto L_x = right.cartesianGTO().center()(CartesianDirection::x);
        const auto L_y = right.cartesianGTO().center()(CartesianDirection::y);
        const auto L_z = right.cartesianGTO().center()(CartesianDirection::z);

        const auto k_K = left.kVector();
        const auto k_L = right.kVector();
        const Vector<double, 3> k1 = right.kVector() - left.kVector();  // The k-vector of the London overlap distribution.

        const auto k_K_x = k_K(CartesianDirection::x);
        const auto k_K_y = k_K(CartesianDirection::y);
        const auto k_K_z = k_K(CartesianDirection::z);

        const auto k_L_x = k_L(CartesianDirection::x);
        const auto k_L_y = k_L(CartesianDirection::y);
        const auto k_L_z = k_L(CartesianDirection::z);

        const auto k1_x = k1(CartesianDirection::x);
        const auto k1_y = k1(CartesianDirection::y);
        const auto k1_z = k1(CartesianDirection::z);


        // For the current component, the integral can be calculated as a product of three contributions.
        PrimitiveOverlapIntegralEngine<LondonGTOShell> S;

        switch (this->component) {
        case CartesianDirection::x: {
            return this->calculate1D(k_K_x, a, K_x, i, k_L_x, b, L_x, j) * S.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case CartesianDirection::y: {
            return S.calculate1D(k1_x, a, K_x, i, b, L_x, j) * this->calculate1D(k_K_y, a, K_y, k, k_L_y, b, L_y, l) * S.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case CartesianDirection::z: {
            return S.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S.calculate1D(k1_y, a, K_y, k, b, L_y, l) * this->calculate1D(k_K_z, a, K_z, m, k_L_z, b, L_z, n);
            break;
        }
        }
    }


    /**
     *  Calculate the linear momentum integral over two London Cartesian GTO 1-D primitives.
     * 
     *  @param k_K              The (directional component of the) k-vector of the left 1-D primitive.
     *  @param a                The Gaussian exponent of the left 1-D primitive.
     *  @param K                The (directional coordinate of the) center of the left 1-D primitive.
     *  @param i                The Cartesian exponent of the left 1-D primitive.
     *  @param k_L              The (directional component of the) k-vector of the right 1-D primitive.
     *  @param b                The Gaussian exponent of the right 1-D primitive.
     *  @param L                The (directional coordinate of the) center of the right 1-D primitive.
     *  @param j                The Cartesian exponent of the right 1-D primitive.
     * 
     *  @return The linear momentum integral over the two London Cartesian GTO given 1-D primitives.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, IntegralScalar> calculate1D(const complex k_K, const double a, const double K, const int i, const complex k_L, const double b, const double L, const int j) {

        using namespace GQCP::literals;

        // The linear momentum integral is a sum of three 1-D overlap integrals. We'll order them from highest to lowest angular momentum.
        const auto k1 = k_L - k_K;  // The (directional component of the) k-vector of the London overlap distribution.
        PrimitiveOverlapIntegralEngine<LondonGTOShell> S;

        return 2.0 * 1.0_ii * b * S.calculate1D(k1, a, K, i, b, L, j + 1) -
               k_L * S.calculate1D(k1, a, K, i, b, L, j) -
               1.0_ii * static_cast<double>(j) * S.calculate1D(k1, a, K, i, b, L, j - 1);
    }
};


}  // namespace GQCP
