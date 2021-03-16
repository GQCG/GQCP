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
#include "Basis/Integrals/Primitive/PrimitiveElectronicDipoleIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveLinearMomentumIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveOverlapIntegralEngine.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Operator/FirstQuantized/AngularMomentumOperator.hpp"


namespace GQCP {


/**
 *  A class that can calculate integrals over the angular momentum operator.
 * 
 *  @tparam _Shell              The type of shell that this integral engine is related to.
 */
template <typename _Shell>
class PrimitiveAngularMomentumIntegralEngine:
    public BaseVectorPrimitiveIntegralEngine {
public:
    // The type of shell that this integral engine is related to.
    using Shell = _Shell;

    // The type of primitive that underlies the type of shell.
    using Primitive = typename Shell::Primitive;

    // The scalar representation of an angular momentum integral.
    using IntegralScalar = product_t<AngularMomentumOperator::Scalar, typename Primitive::Valued>;


private:
    // The angular momentum operator over which the integrals are calculated.
    AngularMomentumOperator angular_momentum_operator;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param angular_momentum_operator                The angular momentum operator over which the integrals are calculated.
     *  @param component                                The initial component of the angular momentum operator over which this engine should calculate integrals.
     */
    PrimitiveAngularMomentumIntegralEngine(const AngularMomentumOperator& angular_momentum_operator, const CartesianDirection component = CartesianDirection::x) :
        angular_momentum_operator {angular_momentum_operator},
        BaseVectorPrimitiveIntegralEngine(component) {}


    /*
    *  MARK: CartesianGTO integrals
    */

    /**
     *  Calculate the angular momentum integral (of the current component) over the two Cartesian GTOs.
     * 
     *  @param left             The left Cartesian GTO.
     *  @param right            The right Cartesian GTO.
     * 
     *  @return The angular momentum integral over the two given Cartesian GTOs.
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


        // For each component of the angular momentum operator, the integrals can be calculated through overlap integrals, linear momentum integrals and position/dipole integrals.
        PrimitiveOverlapIntegralEngine<GTOShell> S0;
        PrimitiveLinearMomentumIntegralEngine<GTOShell> T1;
        PrimitiveElectronicDipoleIntegralEngine<GTOShell> S1 {ElectronicDipoleOperator(this->angular_momentum_operator.reference())};


        // The sign factors in the following formulas are required to switch from electronic dipole integrals to position integrals.
        switch (this->component) {
        case CartesianDirection::x: {
            S1.prepareStateForComponent(CartesianDirection::y);
            T1.prepareStateForComponent(CartesianDirection::z);
            const IntegralScalar term1 = -S1.calculate1D(a, K_y, k, b, L_y, l) * T1.calculate1D(a, K_z, m, b, L_z, n);

            T1.prepareStateForComponent(CartesianDirection::y);
            S1.prepareStateForComponent(CartesianDirection::z);
            const IntegralScalar term2 = -T1.calculate1D(a, K_y, k, b, L_y, l) * S1.calculate1D(a, K_z, m, b, L_z, n);

            return S0.calculate1D(a, K_x, i, b, L_x, j) * (term1 - term2);  // Calculate a component of the cross product.
            break;
        }

        case CartesianDirection::y: {
            S1.prepareStateForComponent(CartesianDirection::z);
            T1.prepareStateForComponent(CartesianDirection::x);
            const IntegralScalar term1 = -S1.calculate1D(a, K_z, m, b, L_z, n) * T1.calculate1D(a, K_x, i, b, L_x, j);

            T1.prepareStateForComponent(CartesianDirection::z);
            S1.prepareStateForComponent(CartesianDirection::x);
            const IntegralScalar term2 = -T1.calculate1D(a, K_z, m, b, L_z, n) * S1.calculate1D(a, K_x, i, b, L_x, j);

            return S0.calculate1D(a, K_y, k, b, L_y, l) * (term1 - term2);  // Calculate a component of the cross product.
            break;
        }

        case CartesianDirection::z: {
            S1.prepareStateForComponent(CartesianDirection::x);
            T1.prepareStateForComponent(CartesianDirection::y);
            const IntegralScalar term1 = -S1.calculate1D(a, K_x, i, b, L_x, j) * T1.calculate1D(a, K_y, k, b, L_y, l);

            T1.prepareStateForComponent(CartesianDirection::x);
            S1.prepareStateForComponent(CartesianDirection::y);
            const IntegralScalar term2 = -T1.calculate1D(a, K_x, i, b, L_x, j) * S1.calculate1D(a, K_y, k, b, L_y, l);

            return S0.calculate1D(a, K_z, m, b, L_z, n) * (term1 - term2);  // Calculate a component of the cross product.
            break;
        }
        }
    }


    /*
     *  MARK: London CartesianGTO integrals
     */

    /**
     *  Calculate the angular momentum integral over two London Cartesian GTOs.
     * 
     *  @param left             The left London Cartesian GTO.
     *  @param right            The right London Cartesian GTO.
     * 
     *  @return The angular momentum integral over the two given London Cartesian GTOs.
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


        // For each component of the angular momentum operator, the integrals can be calculated through overlap integrals, linear momentum integrals and position/dipole integrals.
        PrimitiveOverlapIntegralEngine<LondonGTOShell> S0;
        PrimitiveLinearMomentumIntegralEngine<LondonGTOShell> T1;
        PrimitiveElectronicDipoleIntegralEngine<LondonGTOShell> S1 {ElectronicDipoleOperator(this->angular_momentum_operator.reference())};


        // The sign factors in the following formulas are required to switch from electronic dipole integrals to position integrals.
        switch (this->component) {
        case CartesianDirection::x: {
            S1.prepareStateForComponent(CartesianDirection::y);
            T1.prepareStateForComponent(CartesianDirection::z);
            const IntegralScalar term1 = -S1.calculate1D(k1_y, a, K_y, k, b, L_y, l) * T1.calculate1D(k_K_z, a, K_z, m, k_L_z, b, L_z, n);

            T1.prepareStateForComponent(CartesianDirection::y);
            S1.prepareStateForComponent(CartesianDirection::z);
            const IntegralScalar term2 = -T1.calculate1D(k_K_y, a, K_y, k, k_L_y, b, L_y, l) * S1.calculate1D(k1_z, a, K_z, m, b, L_z, n);

            return S0.calculate1D(k1_x, a, K_x, i, b, L_x, j) * (term1 - term2);  // Calculate a component of the cross product.
            break;
        }

        case CartesianDirection::y: {
            S1.prepareStateForComponent(CartesianDirection::z);
            T1.prepareStateForComponent(CartesianDirection::x);
            const IntegralScalar term1 = -S1.calculate1D(k1_z, a, K_z, m, b, L_z, n) * T1.calculate1D(k_K_x, a, K_x, i, k_L_x, b, L_x, j);

            T1.prepareStateForComponent(CartesianDirection::z);
            S1.prepareStateForComponent(CartesianDirection::x);
            const IntegralScalar term2 = -T1.calculate1D(k_K_z, a, K_z, m, k_L_z, b, L_z, n) * S1.calculate1D(k1_x, a, K_x, i, b, L_x, j);

            return S0.calculate1D(k1_y, a, K_y, k, b, L_y, l) * (term1 - term2);  // Calculate a component of the cross product.
            break;
        }

        case CartesianDirection::z: {
            S1.prepareStateForComponent(CartesianDirection::x);
            T1.prepareStateForComponent(CartesianDirection::y);
            const IntegralScalar term1 = -S1.calculate1D(k1_x, a, K_x, i, b, L_x, j) * T1.calculate1D(k_K_y, a, K_y, k, k_L_y, b, L_y, l);

            T1.prepareStateForComponent(CartesianDirection::x);
            S1.prepareStateForComponent(CartesianDirection::y);
            const IntegralScalar term2 = -T1.calculate1D(k_K_x, a, K_x, i, k_L_x, b, L_x, j) * S1.calculate1D(k1_y, a, K_y, k, b, L_y, l);

            return S0.calculate1D(k1_z, a, K_z, m, b, L_z, n) * (term1 - term2);  // Calculate a component of the cross product.
            break;
        }
        }
    }
};


}  // namespace GQCP
