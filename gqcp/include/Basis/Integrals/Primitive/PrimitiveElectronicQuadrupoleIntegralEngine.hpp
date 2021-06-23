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

#include "Basis/Integrals/Primitive/BaseMatrixPrimitiveIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveElectronicDipoleIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveOverlapIntegralEngine.hpp"
#include "Basis/ScalarBasis/LondonGTOShell.hpp"
#include "Mathematical/Functions/LondonCartesianGTO.hpp"
#include "Operator/FirstQuantized/ElectronicQuadrupoleOperator.hpp"
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
class PrimitiveElectronicQuadrupoleIntegralEngine:
    public BaseMatrixPrimitiveIntegralEngine {
public:
    // The type of shell that this integral engine is related to.
    using Shell = _Shell;

    // The type of primitive that underlies the type of shell.
    using Primitive = typename Shell::Primitive;

    // The scalar representation of an overlap integral.
    using IntegralScalar = product_t<ElectronicQuadrupoleOperator::Scalar, typename Primitive::OutputType>;


private:
    // The electronic quadrupole operator over which integrals should be calculated.
    ElectronicQuadrupoleOperator quadrupole_operator;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param quadrupole_operator          The electronic quadrupole operator over which integrals should be calculated.
     *  @param component                    The initial component of the electronic quadrupole operator over which this primitive engine should calculate integrals.
     */
    PrimitiveElectronicQuadrupoleIntegralEngine(const ElectronicQuadrupoleOperator& quadrupole_operator, const DyadicCartesianDirection component = DyadicCartesianDirection::xx) :
        quadrupole_operator {quadrupole_operator},
        BaseMatrixPrimitiveIntegralEngine(component) {}


    /*
     *  MARK: London CartesianGTO integrals
     */

    /**
     *  Calculate the electronic quadrupole integral over two London Cartesian GTOs.
     * 
     *  @param left             The left London Cartesian GTO.
     *  @param right            The right London Cartesian GTO.
     * 
     *  @return The electronic quadrupole integral over the two given London Cartesian GTOs.
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


        // For the current component, the integral can be calculated as a product of three contributions, containing one-dimensional quadrupole integrals, one-dimensional dipole integrals and one-dimensional overlap integrals.
        // Even though we are working with electronic dipole integrals as opposed to working with position integrals, there is no need to incorporate any sign factors because the dipole integrals always appear in pairs, hence the potential sign factors cancel out.
        PrimitiveOverlapIntegralEngine<LondonGTOShell> S0;
        PrimitiveElectronicDipoleIntegralEngine<LondonGTOShell> S1 {ElectronicDipoleOperator(this->quadrupole_operator.reference())};

        switch (this->component) {
        case DyadicCartesianDirection::xx: {
            return this->calculate1D(k1_x, a, K_x, i, b, L_x, j) * S0.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S0.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case DyadicCartesianDirection::xy: {
            return S1.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S1.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S0.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case DyadicCartesianDirection::xz: {
            return S1.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S0.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S1.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case DyadicCartesianDirection::yx: {
            return S1.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S1.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S0.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case DyadicCartesianDirection::yy: {
            return S0.calculate1D(k1_x, a, K_x, i, b, L_x, j) * this->calculate1D(k1_y, a, K_y, k, b, L_y, l) * S0.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case DyadicCartesianDirection::yz: {
            return S0.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S1.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S1.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case DyadicCartesianDirection::zx: {
            return S1.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S0.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S1.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case DyadicCartesianDirection::zy: {
            return S0.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S1.calculate1D(k1_y, a, K_y, k, b, L_y, l) * S1.calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }

        case DyadicCartesianDirection::zz: {
            return S0.calculate1D(k1_x, a, K_x, i, b, L_x, j) * S0.calculate1D(k1_y, a, K_y, k, b, L_y, l) * this->calculate1D(k1_z, a, K_z, m, b, L_z, n);
            break;
        }
        }
    }


    /**
     *  Calculate the electronic quadrupole integral over two London Cartesian GTO 1-D primitives.
     * 
     *  @param k1               The (directional component of the) k-vector of the London overlap distribution.
     *  @param a                The Gaussian exponent of the left 1-D primitive.
     *  @param K                The (directional coordinate of the) center of the left 1-D primitive.
     *  @param i                The Cartesian exponent of the left 1-D primitive.
     *  @param b                The Gaussian exponent of the right 1-D primitive.
     *  @param L                The (directional coordinate of the) center of the right 1-D primitive.
     *  @param j                The Cartesian exponent of the right 1-D primitive.
     * 
     *  @return The electronic quadrupole integral over the two London Cartesian GTO 1-D primitives.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, IntegralScalar> calculate1D(const complex k1, const double a, const double K, const int i, const double b, const double L, const int j) {

        // Prepare the component X_KC, which is the x-, y- or z-component of the vector between K and the origin of the quadrupole operator.
        CartesianDirection component;
        switch (this->component) {
        case (DyadicCartesianDirection::xx): {
            component = CartesianDirection::x;
            break;
        }
        case (DyadicCartesianDirection::yy): {
            component = CartesianDirection::y;
            break;
        }
        case (DyadicCartesianDirection::zz): {
            component = CartesianDirection::z;
            break;
        }
        default:
            break;
        }
        const auto X_KC = K - this->quadrupole_operator.reference()(component);


        // The 1-D electronic quadrupole integral can be calculated completely from overlap integrals.
        PrimitiveOverlapIntegralEngine<LondonGTOShell> S0;
        return S0.calculate1D(k1, a, K, i + 2, b, L, j) +
               2 * X_KC * S0.calculate1D(k1, a, K, i + 1, b, L, j) +
               std::pow(X_KC, 2) * S0.calculate1D(k1, a, K, i, b, L, j);
    }
};


}  // namespace GQCP
