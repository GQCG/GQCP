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

#include "Basis/Integrals/Primitive/BaseScalarPrimitiveIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/McMurchieDavidsonCoefficient.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/LondonGTOShell.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/LondonCartesianGTO.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Utilities/aliases.hpp"
#include "Utilities/literals.hpp"
#include "Utilities/type_traits.hpp"

#include <boost/math/constants/constants.hpp>


namespace GQCP {


/**
 *  A class that can calculate overlap integrals.
 * 
 *  @tparam _Shell              The type of shell that this integral engine is related to.
 */
template <typename _Shell>
class PrimitiveOverlapIntegralEngine:
    public BaseScalarPrimitiveIntegralEngine {
public:
    // The type of shell that this integral engine is related to.
    using Shell = _Shell;

    // The type of primitive that underlies the type of shell.
    using Primitive = typename Shell::Primitive;

    // The number of components the overlap operator has.
    static constexpr auto Components = OverlapOperator::NumberOfComponents;

    // The scalar representation of an overlap integral.
    using IntegralScalar = product_t<OverlapOperator::Scalar, typename Primitive::Valued>;


public:
    /*
     *  MARK: CartesianGTO integrals
     */

    /**
     *  Calculate the overlap integral over two Cartesian GTOs.
     * 
     *  @param left             The left Cartesian GTO.
     *  @param right            The right Cartesian GTO.
     * 
     *  @return The overlap integral over the two given Cartesian GTOs.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, GTOShell>::value, IntegralScalar> calculate(const CartesianGTO& left, const CartesianGTO& right) {

        // The 3D integral is separable in three 1D integrals.
        IntegralScalar primitive_integral {1.0};
        for (const auto& direction : {GQCP::CartesianDirection::x, GQCP::CartesianDirection::y, GQCP::CartesianDirection::z}) {
            const auto i = left.cartesianExponents().value(direction);
            const auto j = right.cartesianExponents().value(direction);

            primitive_integral *= this->calculate1D(left.gaussianExponent(), left.center()(direction), i, right.gaussianExponent(), right.center()(direction), j);
        }

        return primitive_integral;
    }


    /**
     *  Calculate the overlap integral over two Cartesian GTO 1-D primitives.
     * 
     *  @param a                The Gaussian exponent of the left 1-D primitive.
     *  @param K                The (directional coordinate of the) center of the left 1-D primitive.
     *  @param i                The Cartesian exponent of the left 1-D primitive.
     *  @param b                The Gaussian exponent of the right 1-D primitive.
     *  @param L                The (directional coordinate of the) center of the right 1-D primitive.
     *  @param j                The Cartesian exponent of the right 1-D primitive.
     * 
     *  @return The overlap integral over the two Cartesian GTO 1-D primitives.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, GTOShell>::value, IntegralScalar> calculate1D(const double a, const double K, const int i, const double b, const double L, const int j) {

        // Negative Cartesian exponents should be ignored: the correct value for the corresponding integral is 0.
        if ((i < 0) || (j < 0)) {
            return 0.0;
        }

        // Use the McMurchie-Davidson recursion to calculate the overlap integral.
        const auto p = a + b;
        const McMurchieDavidsonCoefficient E {K, a, L, b};

        return std::pow(boost::math::constants::pi<double>() / p, 0.5) * E(i, j, 0);
    }


    /*
     *  MARK: LondonCartesianGTO integrals
     */

    /**
     *  Calculate the overlap integral over two London Cartesian GTOs.
     * 
     *  @param left             The left London Cartesian GTO.
     *  @param right            The right London Cartesian GTO.
     * 
     *  @return The overlap integral over the two given London Cartesian GTOs.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, IntegralScalar> calculate(const LondonCartesianGTO& left, const LondonCartesianGTO& right) {

        const Vector<double, 3> k1 = right.kVector() - left.kVector();  // The k-vector of the London overlap distribution.

        // The 3D integral is separable in three 1D integrals.
        IntegralScalar primitive_integral {1.0};
        for (const auto& direction : {GQCP::CartesianDirection::x, GQCP::CartesianDirection::y, GQCP::CartesianDirection::z}) {
            const auto a = left.cartesianGTO().gaussianExponent();
            const auto K = left.cartesianGTO().center()(direction);
            const auto i = left.cartesianGTO().cartesianExponents().value(direction);

            const auto b = right.cartesianGTO().gaussianExponent();
            const auto L = right.cartesianGTO().center()(direction);
            const auto j = right.cartesianGTO().cartesianExponents().value(direction);

            const auto k1_component = k1(direction);
            primitive_integral *= this->calculate1D(k1_component, a, K, i, b, L, j);
        }

        return primitive_integral;
    }


    /**
     *  Calculate the overlap integral over two London Cartesian GTO 1-D primitives.
     * 
     *  @param k1               The (directional component of the) k-vector of the London overlap distribution.
     *  @param a                The Gaussian exponent of the left 1-D primitive.
     *  @param K                The (directional coordinate of the) center of the left 1-D primitive.
     *  @param i                The Cartesian exponent of the left 1-D primitive.
     *  @param b                The Gaussian exponent of the right 1-D primitive.
     *  @param L                The (directional coordinate of the) center of the right 1-D primitive.
     *  @param j                The Cartesian exponent of the right 1-D primitive.
     * 
     *  @return The overlap integral over the two London Cartesian GTO 1-D primitives.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, IntegralScalar> calculate1D(const complex k1, const double a, const double K, const int i, const double b, const double L, const int j) {

        // Negative Cartesian exponents should be ignored: the correct value for the corresponding integral is 0.
        if ((i < 0) || (j < 0)) {
            return 0.0;
        }

        using namespace GQCP::literals;


        // Use the McMurchie-Davidson recursion to calculate the overlap integral.
        const auto p = a + b;
        const McMurchieDavidsonCoefficient E {K, a, L, b};
        const auto P = E.centerOfMass();

        IntegralScalar integral {0.0};
        for (int t = 0; t <= i + j; t++) {
            integral += E(i, j, t) *
                        std::pow(-1.0_ii * k1, t) *
                        std::pow(boost::math::constants::pi<double>() / p, 0.5) *
                        std::exp(-1.0_ii * k1 * P) *
                        std::exp(-std::pow(k1, 2) / (4 * p));
        }

        return integral;
    }
};


}  // namespace GQCP
