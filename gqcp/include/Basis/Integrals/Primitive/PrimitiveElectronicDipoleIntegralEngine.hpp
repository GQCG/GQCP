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
#include "Operator/FirstQuantized/ElectronicDipoleOperator.hpp"
#include "Utilities/type_traits.hpp"

#include <boost/math/constants/constants.hpp>


namespace GQCP {


/**
 *  A class that can calculate electronic dipole integrals, i.e. over the negative of the position operator.
 * 
 *  @tparam _Shell              The type of shell that this integral engine is related to.
 * 
 *  @note The integrals that this primitive engine produces include the minus sign due to the charge of the electron.
 */
template <typename _Shell>
class PrimitiveElectronicDipoleIntegralEngine:
    public BaseVectorPrimitiveIntegralEngine {
public:
    // The type of shell that this integral engine is related to.
    using Shell = _Shell;

    // The type of primitive that underlies the type of shell.
    using Primitive = typename Shell::Primitive;

    // The scalar representation of an overlap integral.
    using IntegralScalar = product_t<ElectronicDipoleOperator::Scalar, typename Primitive::OutputType>;


private:
    // The electronic dipole operator over which integrals should be calculated.
    ElectronicDipoleOperator dipole_operator;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param dipole_operator              The electronic dipole operator over which integrals should be calculated.
     *  @param component                    The initial component of the electronic dipole operator over which this primitive engine should calculate integrals.
     */
    PrimitiveElectronicDipoleIntegralEngine(const ElectronicDipoleOperator& dipole_operator, const CartesianDirection component = CartesianDirection::x) :
        dipole_operator {dipole_operator},
        BaseVectorPrimitiveIntegralEngine(component) {}


    /*
     *  MARK: CartesianGTO integrals
     */

    /**
     *  Calculate the electronic dipole integral over two Cartesian GTOs.
     * 
     *  @param left             The left Cartesian GTO.
     *  @param right            The right Cartesian GTO.
     * 
     *  @return The electronic dipole integral over the two given Cartesian GTOs.
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
     *  Calculate the electronic dipole integral over two Cartesian GTO 1-D primitives.
     * 
     *  @param a                The Gaussian exponent of the left 1-D primitive.
     *  @param K                The (directional coordinate of the) center of the left 1-D primitive.
     *  @param i                The Cartesian exponent of the left 1-D primitive.
     *  @param b                The Gaussian exponent of the right 1-D primitive.
     *  @param L                The (directional coordinate of the) center of the right 1-D primitive.
     *  @param j                The Cartesian exponent of the right 1-D primitive.
     * 
     *  @return The electronic dipole integral over the two given 1-D primitives.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, GTOShell>::value, IntegralScalar> calculate1D(const double a, const double K, const int i, const double b, const double L, const int j) {

        // Prepare some variables.
        const McMurchieDavidsonCoefficient E {K, a, L, b};
        const auto P = E.centerOfMass();
        const auto p = a + b;

        const auto X_PC = P - this->dipole_operator.reference()(this->component);  // The distance between P and the origin of the dipole operator.


        // Calculate the dipole integral over the current component. The minus sign comes from the charge of the electron.
        return -std::pow(boost::math::constants::pi<double>() / p, 0.5) * (E(i, j, 1) + X_PC * E(i, j, 0));
    }


    /*
     *  MARK: London CartesianGTO integrals
     */

    /**
     *  Calculate the electronic dipole integral over two London Cartesian GTOs.
     * 
     *  @param left             The left London Cartesian GTO.
     *  @param right            The right London Cartesian GTO.
     * 
     *  @return The electronic dipole integral over the two given London Cartesian GTOs.
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

        const Vector<double, 3> k1 = right.kVector() - left.kVector();  // The k-vector of the London overlap distribution.

        const auto k1_x = k1(CartesianDirection::x);
        const auto k1_y = k1(CartesianDirection::y);
        const auto k1_z = k1(CartesianDirection::z);

        PrimitiveOverlapIntegralEngine<LondonGTOShell> S;


        // For the current component, the integral can be calculated as a product of three contributions.
        switch (this->component) {
        case CartesianDirection::x: {
            return this->calculate1D(k1_x, a, K_x, i, b, L_x, j) * S.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case CartesianDirection::y: {
            return S.calculate1D(k1_x, a, K_x, i, b, L_x, j) * this->calculate1D(k1_y, a, K_y, k, b, L_y, l) * S.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case CartesianDirection::z: {
            return S.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S.calculate1D(k1_y, a, K_y, k, b, L_y, l) * this->calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }
        }
    }


    /**
     *  Calculate the electronic dipole integral over two London Cartesian GTO 1-D primitives.
     * 
     *  @param k1               The (directional component of the) k-vector of the London overlap distribution.
     *  @param a                The Gaussian exponent of the left 1-D primitive.
     *  @param K                The (directional coordinate of the) center of the left 1-D primitive.
     *  @param i                The Cartesian exponent of the left 1-D primitive.
     *  @param b                The Gaussian exponent of the right 1-D primitive.
     *  @param L                The (directional coordinate of the) center of the right 1-D primitive.
     *  @param j                The Cartesian exponent of the right 1-D primitive.
     * 
     *  @return The electronic dipole integral over the two London Cartesian GTO 1-D primitives.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, IntegralScalar> calculate1D(const complex k1, const double a, const double K, const int i, const double b, const double L, const int j) {

        // Prepare some variables.
        const auto X_KC = K - this->dipole_operator.reference()(this->component);  // The distance between K and the origin of the dipole operator.


        // The 1-D electronic dipole integral can be calculated completely from overlap integrals. The sign factor is included to account for the sign of the electron.
        PrimitiveOverlapIntegralEngine<LondonGTOShell> S;
        return (-1.0) * (S.calculate1D(k1, a, K, i + 1, b, L, j) +
                         X_KC * S.calculate1D(k1, a, K, i, b, L, j));
    }
};


}  // namespace GQCP
