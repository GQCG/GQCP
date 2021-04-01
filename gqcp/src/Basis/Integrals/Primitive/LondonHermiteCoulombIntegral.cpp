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

#include "Basis/Integrals/Primitive/LondonHermiteCoulombIntegral.hpp"

#include "Mathematical/Functions/BoysFunction.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param k1           The k-vector of the London overlap distribution.
 *  @param p            The exponent of the Hermite Gaussian.
 *  @param P            The center of the Hermite Gaussian.
 *  @param C            The center of the Coulomb potential.
 */
LondonHermiteCoulombIntegral::LondonHermiteCoulombIntegral(const Vector<double, 3>& k1, const double p, const Vector<double, 3>& P, const Vector<double, 3>& C) :
    p {p},
    k1 {k1},
    P {P},
    C {C} {}


/*
 *  MARK: London Hermite Coulomb integral implementation
 */

/**
 *  Calculate the value for the (auxiliary) London Hermite Coulomb integral R^{k1, n}_{tuv}(p, P, C).
 * 
 *  @param n            The order of the London Hermite Coulomb integral, i.e. the order of the Boys function.
 *  @param t            The derivative degree in P_x.
 *  @param u            The derivative degree in P_y.
 *  @param v            The derivative degree in P_z.
 */
complex LondonHermiteCoulombIntegral::operator()(const size_t n, const int t, const int u, const int v) const {

    using namespace GQCP::literals;

    // Prepare some variables.
    const Vector<complex, 3> P_ = this->P - (1.0_ii / (2 * this->p)) * this->k1;  // The modified overlap distribution center.
    const Vector<complex, 3> R_P_C = P_ - this->C;


    // If any of the degrees is smaller than 0, the London Hermite Coulomb integral should vanish.
    if ((t < 0) || (u < 0) || (v < 0)) {
        return 0.0;
    }


    // Provide the base case for (t == u == v == 0).
    if ((t == 0) && (u == 0) && (v == 0)) {
        const complex R2_P_C = R_P_C.array().square().sum();

        return std::pow(-2.0 * this->p, n) *
               std::exp(-1.0_ii * this->k1.dot(this->P)) *
               BoysFunction()(n, p * R2_P_C);
    }


    // Recurrence for v.
    else if ((t == 0) && (u == 0)) {
        return -1.0_ii * this->k1(CartesianDirection::z) * this->operator()(n, t, u, v - 1) +
               static_cast<double>(v - 1) * this->operator()(n + 1, t, u, v - 2) +
               R_P_C(CartesianDirection::z) * this->operator()(n + 1, t, u, v - 1);
    }


    // Recurrence for u.
    else if (t == 0) {
        return -1.0_ii * this->k1(CartesianDirection::y) * this->operator()(n, t, u - 1, v) +
               static_cast<double>(u - 1) * this->operator()(n + 1, t, u - 2, v) +
               R_P_C(CartesianDirection::y) * this->operator()(n + 1, t, u - 1, v);
    }

    // Recurrence for t.
    else {
        return -1.0_ii * this->k1(CartesianDirection::x) * this->operator()(n, t - 1, u, v) +
               static_cast<double>(t - 1) * this->operator()(n + 1, t - 2, u, v) +
               R_P_C(CartesianDirection::x) * this->operator()(n + 1, t - 1, u, v);
    }
}


}  // namespace GQCP
