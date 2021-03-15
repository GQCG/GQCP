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

#include "Basis/Integrals/Primitive/HermiteCoulombIntegral.hpp"
#include "Basis/Integrals/Primitive/McMurchieDavidsonCoefficient.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Operator/FirstQuantized/NuclearAttractionOperator.hpp"

#include <boost/math/constants/constants.hpp>


namespace GQCP {


/**
 *  A class that can calculate nuclear attraction integrals.
 * 
 *  @tparam _Shell              The type of shell that this integral engine is related to.
 */
template <typename _Shell>
class PrimitiveNuclearAttractionIntegralEngine {
public:
    // The type of shell that this integral engine is related to.
    using Shell = _Shell;

    // The type of primitive that underlies the type of shell.
    using Primitive = typename Shell::Primitive;

    // The number of components the nuclear attraction operator has.
    static constexpr auto Components = NuclearAttractionOperator::NumberOfComponents;

    // The scalar representation of a nuclear attraction integral.
    using IntegralScalar = product_t<NuclearAttractionOperator::Scalar, typename Primitive::Valued>;


private:
    // The nuclear attraction operator that this engine can calculate integrals over.
    NuclearAttractionOperator nuclear_attraction_operator;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param nuclear_attraction_operator              The nuclear attraction operator that this engine can calculate integrals over.
     */
    PrimitiveNuclearAttractionIntegralEngine(const NuclearAttractionOperator& nuclear_attraction_operator) :
        nuclear_attraction_operator {nuclear_attraction_operator} {}


    /**
     *  MARK: Components
     */

    /**
     *  Prepare this engine's internal state such that it is able to calculate integrals over the given component of the operator.
     * 
     *  @param component                The index of the component of the operator.
     * 
     *  @note Since the nuclear attraction operator has only 1 component, this method has no effect.
     */
    void prepareStateForComponent(const size_t component) {};


    /*
     *  MARK: CartesianGTO integrals
     */

    /**
     *  Calculate the nuclear attraction integral over two Cartesian GTOs.
     * 
     *  @param left             The left Cartesian GTO.
     *  @param right            The right Cartesian GTO.
     * 
     *  @return The nuclear integral over the two given Cartesian GTOs.
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


        // Prepare the McMurchie-Davidson coefficients.
        const McMurchieDavidsonCoefficient E_x {K_x, a, L_x, b};
        const McMurchieDavidsonCoefficient E_y {K_y, a, L_y, b};
        const McMurchieDavidsonCoefficient E_z {K_z, a, L_z, b};

        const double p = a + b;
        const Vector<double, 3> P {E_x.centerOfMass(), E_y.centerOfMass(), E_z.centerOfMass()};


        // Calculate the contributions from every nuclear center.
        double total_integral {0.0};

        const auto& nuclei = this->nuclear_attraction_operator.nuclearFramework().nucleiAsVector();
        for (const auto& nucleus : nuclei) {
            double integral {0.0};

            const auto& C = nucleus.position();
            const HermiteCoulombIntegral R {p, P, C};

            for (int t = 0; t <= i + j; t++) {
                for (int u = 0; u <= k + l; u++) {
                    for (int v = 0; v <= m + n; v++) {
                        // Add the contribution to the integral. The prefactor will be applied at the end.
                        integral += E_x(i, j, t) * E_y(k, l, u) * E_z(m, n, v) * R(0, t, u, v);
                    }
                }
            }
            const auto charge = static_cast<double>(nucleus.charge());
            total_integral += (-charge) * integral;
        }

        return 2 * boost::math::constants::pi<double>() / p * total_integral;
    }
};


}  // namespace GQCP
