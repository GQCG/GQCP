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

#include "Basis/Integrals/HermiteCoulombIntegral.hpp"

#include "Mathematical/Functions/BoysFunction.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param p            The exponent of the Hermite Gaussian.
 *  @param P            The center of the Hermite Gaussian.
 *  @param C            The nuclear center.
 */
HermiteCoulombIntegral::HermiteCoulombIntegral(const double p, const Vector<double, 3>& P, const Vector<double, 3>& C) :
    p {p},
    P {P},
    C {C} {}


/*
 *  MARK: Hermite Coulomb integral implementation
 */

/**
 *  Calculate the value for the (auxiliary) Hermite Coulomb integral R^n_{tuv}(p, P, C).
 * 
 *  @param n            The order of the Hermite Coulomb integral, i.e. the order of the Boys function.
 *  @param t            The degree of the Hermite polynomial in x.
 *  @param u            The degree of the Hermite polynomial in y.
 *  @param v            The degree of the Hermite polynomial in z.
 */
double HermiteCoulombIntegral::operator()(const size_t n, const int t, const int u, const int v) const {

    const Vector<double, 3> R_PC = this->P - this->C;


    // If any of the arguments is smaller than 0, the Hermite Coulomb integral should vanish.
    if ((t < 0) || (u < 0) || (v < 0)) {
        return 0.0;
    }


    // Provide the base case for (t == u == v == 0).
    if ((t == 0) && (u == 0) && (v == 0)) {
        const double R2_PC = R_PC.squaredNorm();

        return std::pow(-2.0 * this->p, n) * BoysFunction()(n, p * R2_PC);
    }


    // Recurrence for v.
    else if ((t == 0) && (u == 0)) {
        return (v - 1) * this->operator()(n + 1, t, u, v - 2) +
               R_PC(CartesianDirection::z) * this->operator()(n + 1, t, u, v - 1);
    }


    // Recurrence for u.
    else if (t == 0) {
        return (u - 1) * this->operator()(n + 1, t, u - 2, v) +
               R_PC(CartesianDirection::y) * this->operator()(n + 1, t, u - 1, v);
    }

    // Recurrence for t.
    else {
        return (t - 1) * this->operator()(n + 1, t - 2, u, v) +
               R_PC(CartesianDirection::x) * this->operator()(n + 1, t - 1, u, v);
    }
}


}  // namespace GQCP
