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
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
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
     *  @note Since the kinetic energy operator has only 1 component, this method has no effect.
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


        // The 3D kinetic energy integral is a sum of three contributions (dx^2, dy^2, dz^2).
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

        // The kinetic 1D integral is a sum of three 1D overlap integrals.
        PrimitiveOverlapIntegralEngine<GTOShell> primitive_overlap_engine;

        return -2 * std::pow(b, 2) * primitive_overlap_engine.calculate1D(a, K, i, b, L, j + 2) +
               b * (2 * j + 1) * primitive_overlap_engine.calculate1D(a, K, i, b, L, j) -
               0.5 * j * (j - 1) * primitive_overlap_engine.calculate1D(a, K, i, b, L, j - 2);
    }


    /*
     *  MARK: LondonCartesianGTO integrals
     */
};


}  // namespace GQCP
