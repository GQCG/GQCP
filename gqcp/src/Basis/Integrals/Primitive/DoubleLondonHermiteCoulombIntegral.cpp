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

#include "Basis/Integrals/Primitive/DoubleLondonHermiteCoulombIntegral.hpp"

#include "Mathematical/Functions/BoysFunction.hpp"
#include "Utilities/literals.hpp"

namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param k1           The k-vector of the left London overlap distribution.
 *  @param p            The exponent of the left London overlap distribution.
 *  @param P            The center of the left London overlap distribution.
 *  @param k2           The k-vector of the right London overlap distribution.
 *  @param q            The exponent of the right London overlap distribution.
 *  @param Q            The center of the right London overlap distribution.
 */
DoubleLondonHermiteCoulombIntegral::DoubleLondonHermiteCoulombIntegral(const Vector<double, 3>& k1, const double p, const Vector<double, 3>& P, const Vector<double, 3>& k2, const double q, const Vector<double, 3>& Q) :
    p {p},
    q {q},
    k1 {k1},
    k2 {k2},
    P {P},
    Q {Q} {}


/*
 *  MARK: London Hermite Coulomb integral implementation
 */

/**
 *  Calculate the value for the double (auxiliary) London Hermite Coulomb integral R^{k1, k2, n}_{tuv, tau mu nu}(p, q, P, Q).
 * 
 *  @param n            The order of the London Hermite Coulomb integral, i.e. the order of the Boys function.
 *  @param t            The derivative degree in P_x.
 *  @param u            The derivative degree in P_y.
 *  @param v            The derivative degree in P_z.
 *  @param tau          The derivative degree in Q_x.
 *  @param mu           The derivative degree in Q_y.
 *  @param nu           The derivative degree in Q_z.
 */
complex DoubleLondonHermiteCoulombIntegral::operator()(const size_t n, const int t, const int u, const int v, const int tau, const int mu, const int nu) const {

    using namespace GQCP::literals;


    // Prepare some variables.
    const Vector<complex, 3> P_ = this->P - (1.0_ii / (2 * this->p)) * this->k1;  // The left modified overlap distribution center.
    const Vector<complex, 3> Q_ = this->Q - (1.0_ii / (2 * this->q)) * this->k2;  // The right modified overlap distribution center.
    const Vector<complex, 3> R_P_Q_ = P_ - Q_;


    // If any of the degrees is smaller than 0, the London Hermite Coulomb integral should vanish.
    if ((t < 0) || (u < 0) || (v < 0) || (tau < 0) || (mu < 0) || (nu < 0)) {
        return 0.0;
    }


    // Provide the base case for (t == u == v == tau == mu == nu == 0).
    if ((t == 0) && (u == 0) && (v == 0) && (tau == 0) && (mu == 0) && (nu == 0)) {
        const auto alpha = this->p * this->q / (this->p + this->q);
        const complex R2_P_Q_ = R_P_Q_.array().square().sum();

        return std::pow(-2.0 * alpha, n) *
               std::exp(-1.0_ii * this->k1.dot(this->P) - 1.0_ii * this->k2.dot(this->Q)) *
               BoysFunction()(n, alpha * R2_P_Q_);
    }


    // Recurrence for nu.
    else if ((t == 0) && (u == 0) && (v == 0) && (tau == 0) && (mu == 0)) {
        return -1.0_ii * this->k2(CartesianDirection::z) * this->operator()(n, t, u, v, tau, mu, nu - 1) -
               R_P_Q_(CartesianDirection::z) * this->operator()(n + 1, t, u, v, tau, mu, nu - 1) -
               static_cast<double>(v) * this->operator()(n + 1, t, u, v - 1, tau, mu, nu - 1) +
               static_cast<double>(nu - 1) * this->operator()(n + 1, t, u, v, tau, mu, nu - 2);
    }


    // Recurrence for mu.
    else if ((t == 0) && (u == 0) && (v == 0) && (tau == 0)) {
        return -1.0_ii * this->k2(CartesianDirection::y) * this->operator()(n, t, u, v, tau, mu - 1, nu) -
               R_P_Q_(CartesianDirection::y) * this->operator()(n + 1, t, u, v, tau, mu - 1, nu) -
               static_cast<double>(u) * this->operator()(n + 1, t, u - 1, v, tau, mu - 1, nu) +
               static_cast<double>(mu - 1) * this->operator()(n + 1, t, u, v, tau, mu - 2, nu);
    }


    // Recurrence for tau.
    else if ((t == 0) && (u == 0) && (v == 0)) {
        return -1.0_ii * this->k2(CartesianDirection::x) * this->operator()(n, t, u, v, tau - 1, mu, nu) -
               R_P_Q_(CartesianDirection::x) * this->operator()(n + 1, t, u, v, tau - 1, mu, nu) -
               static_cast<double>(t) * this->operator()(n + 1, t - 1, u, v, tau - 1, mu, nu) +
               static_cast<double>(tau - 1) * this->operator()(n + 1, t, u, v, tau - 2, mu, nu);
    }


    // Recurrence for v.
    else if ((t == 0) && (u == 0)) {
        return -1.0_ii * this->k1(CartesianDirection::z) * this->operator()(n, t, u, v - 1, tau, mu, nu) +
               R_P_Q_(CartesianDirection::z) * this->operator()(n + 1, t, u, v - 1, tau, mu, nu) +
               static_cast<double>(v - 1) * this->operator()(n + 1, t, u, v - 2, tau, mu, nu) -
               static_cast<double>(nu) * this->operator()(n + 1, t, u, v - 1, tau, mu, nu - 1);
    }


    // Recurrence for u.
    else if (t == 0) {
        return -1.0_ii * this->k1(CartesianDirection::y) * this->operator()(n, t, u - 1, v, tau, mu, nu) +
               R_P_Q_(CartesianDirection::y) * this->operator()(n + 1, t, u - 1, v, tau, mu, nu) +
               static_cast<double>(u - 1) * this->operator()(n + 1, t, u - 2, v, tau, mu, nu) -
               static_cast<double>(mu) * this->operator()(n + 1, t, u - 1, v, tau, mu - 1, nu);
    }


    // Recurrence for t.
    else {
        return -1.0_ii * this->k1(CartesianDirection::x) * this->operator()(n, t - 1, u, v, tau, mu, nu) +
               R_P_Q_(CartesianDirection::x) * this->operator()(n + 1, t - 1, u, v, tau, mu, nu) +
               static_cast<double>(t - 1) * this->operator()(n + 1, t - 2, u, v, tau, mu, nu) -
               static_cast<double>(tau) * this->operator()(n + 1, t - 1, u, v, tau - 1, mu, nu);
    }
}


}  // namespace GQCP
