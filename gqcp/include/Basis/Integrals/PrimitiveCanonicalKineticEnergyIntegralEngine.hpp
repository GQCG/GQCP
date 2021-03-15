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

#include "Basis/Integrals/PrimitiveOverlapIntegralEngine.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/LondonGTOShell.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/LondonCartesianGTO.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
#include "Utilities/literals.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A class that can calculate integrals over the canonical kinetic energy operator, i.e. - 1/2 nabla^2.
 * 
 *  @tparam _Shell              The type of shell that this integral engine is related to.
 */
template <typename _Shell>
class PrimitiveCanonicalKineticEnergyIntegralEngine {
public:
    // The type of shell that this integral engine is related to.
    using Shell = _Shell;

    // The type of primitive that underlies the type of shell.
    using Primitive = typename Shell::Primitive;

    // The number of components the canonical kinetic energy operator has.
    static constexpr auto Components = KineticOperator::NumberOfComponents;

    // The scalar representation of a canonical kinetic energy integral.
    using IntegralScalar = product_t<KineticOperator::Scalar, typename Primitive::Valued>;


public:
    /**
     *  MARK: Components
     */

    /**
     *  Prepare this engine's internal state such that it is able to calculate integrals over the given component of the operator.
     * 
     *  @param component                the index of the component of the operator
     * 
     *  @note Since the canonical kinetic energy operator has only 1 component, this method has no effect.
     */
    void prepareStateForComponent(const size_t component) {};


    /*
     *  MARK: CartesianGTO integrals
     */

    /**
     *  Calculate the canonical kinetic energy integral over two Cartesian GTOs.
     * 
     *  @param left             The left Cartesian GTO.
     *  @param right            The right Cartesian GTO.
     * 
     *  @return The canonical kinetic energy integral over the two given Cartesian GTOs.
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


        // The 3D canonical kinetic energy integral is a sum of three contributions (dx^2, dy^2, dz^2).
        PrimitiveOverlapIntegralEngine<GTOShell> primitive_overlap_engine;

        IntegralScalar primitive_integral = 1.0;
        return this->calculate1D(a, K_x, i, b, L_x, j) * primitive_overlap_engine.calculate1D(a, K_y, k, b, L_y, l) * primitive_overlap_engine.calculate1D(a, K_z, m, b, L_z, n) +
               primitive_overlap_engine.calculate1D(a, K_x, i, b, L_x, j) * this->calculate1D(a, K_y, k, b, L_y, l) * primitive_overlap_engine.calculate1D(a, K_z, m, b, L_z, n) +
               primitive_overlap_engine.calculate1D(a, K_x, i, b, L_x, j) * primitive_overlap_engine.calculate1D(a, K_y, k, b, L_y, l) * this->calculate1D(a, K_z, m, b, L_z, n);
    }


    /**
     *  Calculate the canonical kinetic energy integral over two Cartesian GTO 1-D primitives.
     * 
     *  @param a                The Gaussian exponent of the left 1-D primitive.
     *  @param K                The (directional coordinate of the) center of the left 1-D primitive.
     *  @param i                The Cartesian exponent of the left 1-D primitive.
     *  @param b                The Gaussian exponent of the right 1-D primitive.
     *  @param L                The (directional coordinate of the) center of the right 1-D primitive.
     *  @param j                The Cartesian exponent of the right 1-D primitive.
     * 
     *  @return The canonical kinetic energy integral over the two Cartesian GTO given 1-D primitives.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, GTOShell>::value, IntegralScalar> calculate1D(const double a, const double K, const int i, const double b, const double L, const int j) {

        // The canonical kinetic 1D integral is a sum of three 1D overlap integrals.
        PrimitiveOverlapIntegralEngine<GTOShell> primitive_overlap_engine;

        return -2 * std::pow(b, 2) * primitive_overlap_engine.calculate1D(a, K, i, b, L, j + 2) +
               b * (2 * j + 1) * primitive_overlap_engine.calculate1D(a, K, i, b, L, j) -
               0.5 * j * (j - 1) * primitive_overlap_engine.calculate1D(a, K, i, b, L, j - 2);
    }


    /*
     *  MARK: LondonCartesianGTO integrals
     */

    /**
     *  Calculate the canonical kinetic energy integral over two London Cartesian GTOs.
     * 
     *  @param left             The left London Cartesian GTO.
     *  @param right            The right London Cartesian GTO.
     * 
     *  @return The canonical kinetic energy integral over the two given London Cartesian GTOs.
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


        // The 3D canonical kinetic energy integral is a sum of three contributions (dx^2, dy^2, dz^2).
        PrimitiveOverlapIntegralEngine<LondonGTOShell> S;

        IntegralScalar primitive_integral = 1.0;
        return this->calculate1D(k_K_x, a, K_x, i, k_L_x, b, L_x, j) * S.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S.calculate1D(k1_z, a, K_z, m, b, L_z, n) +
               S.calculate1D(k1_x, a, K_x, i, b, L_x, j) * this->calculate1D(k_K_y, a, K_y, k, k_L_y, b, L_y, l) * S.calculate1D(k1_z, a, K_z, m, b, L_z, n) +
               S.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S.calculate1D(k1_y, a, K_y, k, b, L_y, l) * this->calculate1D(k_K_z, a, K_z, m, k_L_z, b, L_z, n);
    }


    /**
     *  Calculate the canonical kinetic energy integral over two London Cartesian GTO 1-D primitives.
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
     *  @return The canonical kinetic energy integral over the two London Cartesian GTO given 1-D primitives.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, IntegralScalar> calculate1D(const complex k_K, const double a, const double K, const int i, const complex k_L, const double b, const double L, const int j) {

        using namespace GQCP::literals;

        // The canonical kinetic 1D integral is a sum of five 1-D overlap integrals. We'll order them from highest to lowest angular momentum.
        const auto k1 = k_L - k_K;  // The (directional component of the) k-vector of the London overlap distribution.
        PrimitiveOverlapIntegralEngine<LondonGTOShell> S;

        return -2 * std::pow(b, 2) * S.calculate1D(k1, a, K, i, b, L, j + 2) -
               2 * b * 1.0_ii * k_L * S.calculate1D(k1, a, K, i, b, L, j + 1) +
               (b * (2 * j + 1) + 0.5 * std::pow(k_L, 2)) * S.calculate1D(k1, a, K, i, b, L, j) +
               static_cast<double>(j) * 1.0_ii * k_L * S.calculate1D(k1, a, K, i, b, L, j - 1) -
               0.5 * j * (j - 1) * S.calculate1D(k1, a, K, i, b, L, j - 2);
    }
};


}  // namespace GQCP
