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

#include "Basis/Integrals/PrimitiveKineticEnergyIntegralEngine.hpp"

#include "Basis/Integrals/PrimitiveOverlapIntegralEngine.hpp"


namespace GQCP {


/**
 *  @param left             the left Cartesian GTO (primitive)
 *  @param right            the right Cartesian GTO (primitive)
 * 
 *  @return the kinetic energy integral over the two given primitives
 */
double PrimitiveKineticEnergyIntegralEngine::calculate(const CartesianGTO& left, const CartesianGTO& right) {

    // Prepare some variables.
    const auto i = static_cast<int>(left.cartesianExponents().value(CartesianDirection::x));
    const auto k = static_cast<int>(left.cartesianExponents().value(CartesianDirection::y));
    const auto m = static_cast<int>(left.cartesianExponents().value(CartesianDirection::z));

    const auto j = static_cast<int>(right.cartesianExponents().value(CartesianDirection::x));
    const auto l = static_cast<int>(right.cartesianExponents().value(CartesianDirection::y));
    const auto n = static_cast<int>(right.cartesianExponents().value(CartesianDirection::z));

    const auto alpha = left.gaussianExponent();
    const auto beta = right.gaussianExponent();

    const auto K_x = left.center()(CartesianDirection::x);
    const auto K_y = left.center()(CartesianDirection::y);
    const auto K_z = left.center()(CartesianDirection::z);

    const auto L_x = right.center()(CartesianDirection::x);
    const auto L_y = right.center()(CartesianDirection::y);
    const auto L_z = right.center()(CartesianDirection::z);


    // The 3D kinetic energy integral is a sum of three contributions (dx^2, dy^2, dz^2).
    PrimitiveOverlapIntegralEngine primitive_overlap_engine;

    double primitive_integral = 1.0;
    return this->calculate1D(alpha, K_x, i, beta, L_x, j) * primitive_overlap_engine.calculate1D(alpha, K_y, k, beta, L_y, l) * primitive_overlap_engine.calculate1D(alpha, K_z, m, beta, L_z, n) +
           primitive_overlap_engine.calculate1D(alpha, K_x, i, beta, L_x, j) * this->calculate1D(alpha, K_y, k, beta, L_y, l) * primitive_overlap_engine.calculate1D(alpha, K_z, m, beta, L_z, n) +
           primitive_overlap_engine.calculate1D(alpha, K_x, i, beta, L_x, j) * primitive_overlap_engine.calculate1D(alpha, K_y, k, beta, L_y, l) * this->calculate1D(alpha, K_z, m, beta, L_z, n);
}

/**
 *  @param alpha            the Gaussian exponent of the left 1-D primitive
 *  @param K                the (directional coordinate of the) center of the left 1-D primitive
 *  @param i                the Cartesian exponent of the left 1-D primitive
 *  @param beta             the Gaussian exponent of the right 1-D primitive
 *  @param L                the (directional coordinate of the) center of the right 1-D primitive
 *  @param j                the Cartesian exponent of the right 1-D primitive
 * 
 *  @return the kinetic energy integral over the two given 1-D primitives
 */
double PrimitiveKineticEnergyIntegralEngine::calculate1D(const double alpha, const double K, const int i, const double beta, const double L, const int j) {

    // The kinetic 1D integral is a sum of three 1D overlap integrals.
    PrimitiveOverlapIntegralEngine primitive_overlap_engine;

    return -2 * std::pow(beta, 2) * primitive_overlap_engine.calculate1D(alpha, K, i, beta, L, j + 2) +
           beta * (2 * j + 1) * primitive_overlap_engine.calculate1D(alpha, K, i, beta, L, j) -
           0.5 * j * (j - 1) * primitive_overlap_engine.calculate1D(alpha, K, i, beta, L, j - 2);
}


}  // namespace GQCP
