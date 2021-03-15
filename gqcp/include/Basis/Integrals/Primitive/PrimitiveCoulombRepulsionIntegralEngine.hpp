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
#include "Basis/Integrals/Primitive/HermiteCoulombIntegral.hpp"
#include "Basis/Integrals/Primitive/McMurchieDavidsonCoefficient.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"

#include <boost/math/constants/constants.hpp>

#include <cmath>


namespace GQCP {


/**
 *  A class that can calculate electron repulsion integrals.
 * 
 *  @tparam _Shell              The type of shell that this integral engine is related to.
 */
template <typename _Shell>
class PrimitiveCoulombRepulsionIntegralEngine:
    public BaseScalarPrimitiveIntegralEngine {
public:
    // The type of shell that this integral engine is related to.
    using Shell = _Shell;

    // The type of primitive that underlies the type of shell.
    using Primitive = typename Shell::Primitive;

    // The scalar representation of a nuclear attraction integral.
    using IntegralScalar = product_t<CoulombRepulsionOperator::Scalar, typename Primitive::Valued>;


public:
    /*
     *  MARK: CartesianGTO integrals
     */

    /**
     *  Calculate the Coulomb repulsion integral over four Cartesian GTOs.
     * 
     *  @param left1            The first left Cartesian GTO.
     *  @param left2            The second left Cartesian GTO.
     *  @param right1           The first right Cartesian GTO.
     *  @param right2           The second right Cartesian GTO.
     * 
     *  @return The Coulomb repulsion integral over the four given Cartesian GTOs.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, GTOShell>::value, IntegralScalar> calculate(const CartesianGTO& left1, const CartesianGTO& left2, const CartesianGTO& right1, const CartesianGTO& right2) {

        // Prepare some variables. Those with an extra underscore represent the 'primed' indices in the notes.
        const auto i = static_cast<int>(left1.cartesianExponents().value(CartesianDirection::x));
        const auto k = static_cast<int>(left1.cartesianExponents().value(CartesianDirection::y));
        const auto m = static_cast<int>(left1.cartesianExponents().value(CartesianDirection::z));

        const auto j = static_cast<int>(left2.cartesianExponents().value(CartesianDirection::x));
        const auto l = static_cast<int>(left2.cartesianExponents().value(CartesianDirection::y));
        const auto n = static_cast<int>(left2.cartesianExponents().value(CartesianDirection::z));

        const auto i_ = static_cast<int>(right1.cartesianExponents().value(CartesianDirection::x));
        const auto k_ = static_cast<int>(right1.cartesianExponents().value(CartesianDirection::y));
        const auto m_ = static_cast<int>(right1.cartesianExponents().value(CartesianDirection::z));

        const auto j_ = static_cast<int>(right2.cartesianExponents().value(CartesianDirection::x));
        const auto l_ = static_cast<int>(right2.cartesianExponents().value(CartesianDirection::y));
        const auto n_ = static_cast<int>(right2.cartesianExponents().value(CartesianDirection::z));

        const auto a = left1.gaussianExponent();
        const auto b = left2.gaussianExponent();
        const auto c = right1.gaussianExponent();
        const auto d = right2.gaussianExponent();

        const auto K_x = left1.center()(CartesianDirection::x);
        const auto K_y = left1.center()(CartesianDirection::y);
        const auto K_z = left1.center()(CartesianDirection::z);

        const auto L_x = left2.center()(CartesianDirection::x);
        const auto L_y = left2.center()(CartesianDirection::y);
        const auto L_z = left2.center()(CartesianDirection::z);

        const auto M_x = right1.center()(CartesianDirection::x);
        const auto M_y = right1.center()(CartesianDirection::y);
        const auto M_z = right1.center()(CartesianDirection::z);

        const auto N_x = right2.center()(CartesianDirection::x);
        const auto N_y = right2.center()(CartesianDirection::y);
        const auto N_z = right2.center()(CartesianDirection::z);


        // Prepare the McMurchie-Davidson coefficients.
        const McMurchieDavidsonCoefficient E_x {K_x, a, L_x, b};
        const McMurchieDavidsonCoefficient E_y {K_y, a, L_y, b};
        const McMurchieDavidsonCoefficient E_z {K_z, a, L_z, b};

        const McMurchieDavidsonCoefficient E_x_ {M_x, c, N_x, d};
        const McMurchieDavidsonCoefficient E_y_ {M_y, c, N_y, d};
        const McMurchieDavidsonCoefficient E_z_ {M_z, c, N_z, d};


        // Prepare the Hermite Coulomb integral.
        const double p = a + b;
        const double q = c + d;
        const double alpha = p * q / (p + q);

        const Vector<double, 3> P {E_x.centerOfMass(), E_y.centerOfMass(), E_z.centerOfMass()};
        const Vector<double, 3> Q {E_x_.centerOfMass(), E_y_.centerOfMass(), E_z_.centerOfMass()};

        const HermiteCoulombIntegral R {alpha, P, Q};


        // Calculate the Coulomb repulsion integrals over the primitives.
        double integral {};
        for (int t = 0; t <= i + j; t++) {
            for (int u = 0; u <= k + l; u++) {
                for (int v = 0; v <= m + n; v++) {
                    for (int tau = 0; tau <= i_ + j_; tau++) {
                        for (int mu = 0; mu <= k_ + l_; mu++) {
                            for (int nu = 0; nu <= m_ + n_; nu++) {
                                // Add the contribution to the integral. The prefactor will be applied at the end.
                                integral += E_x(i, j, t) * E_y(k, l, u) * E_z(m, n, v) *
                                            E_x_(i_, j_, tau) * E_y_(k_, l_, mu) * E_z_(m_, n_, nu) *
                                            std::pow(-1, tau + mu + nu) * R(0, t + tau, u + mu, v + nu);
                            }
                        }
                    }
                }
            }
        }

        return 2 * std::pow(boost::math::constants::pi<double>(), 2.5) / (p * q * std::sqrt(p + q)) * integral;
    }
};


}  // namespace GQCP
