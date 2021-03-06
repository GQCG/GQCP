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

#include "Basis/Integrals/PrimitiveLinearMomentumIntegralEngine.hpp"

#include "Basis/Integrals/PrimitiveOverlapIntegralEngine.hpp"
#include "Utilities/literals.hpp"

namespace GQCP {


/*
 *  PUBLIC METHODS
 */

/**
 *  @param left             the left Cartesian GTO (primitive)
 *  @param right            the right Cartesian GTO (primitive)
 * 
 *  @return the linear momentum integral (of the current component) over the two given primitives
 */
PrimitiveLinearMomentumIntegralEngine::IntegralScalar PrimitiveLinearMomentumIntegralEngine::calculate(const CartesianGTO& left, const CartesianGTO& right) {

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

    PrimitiveOverlapIntegralEngine overlap_engine;


    // For the current component, the integral can be calculated as a product of three contributions.
    switch (this->component) {
    case CartesianDirection::x: {
        return this->calculate1D(alpha, K_x, i, beta, L_x, j) * overlap_engine.calculate1D(alpha, K_y, k, beta, L_y, l) * overlap_engine.calculate1D(alpha, K_z, m, beta, L_z, n);
        break;
    }

    case CartesianDirection::y: {
        return overlap_engine.calculate1D(alpha, K_x, i, beta, L_x, j) * this->calculate1D(alpha, K_y, k, beta, L_y, l) * overlap_engine.calculate1D(alpha, K_z, m, beta, L_z, n);
        break;
    }

    case CartesianDirection::z: {
        return overlap_engine.calculate1D(alpha, K_x, i, beta, L_x, j) * overlap_engine.calculate1D(alpha, K_y, k, beta, L_y, l) * this->calculate1D(alpha, K_z, m, beta, L_z, n);
        break;
    }
    }
}


/**
 *  @param alpha            the Gaussian exponent of the left 1-D primitive
 *  @param K                the (directional coordinate of the) center of the left 1-D primitive
 *  @param i                the Cartesian exponent of the left 1-D primitive
 *  @param beta             the Gaussian exponent of the right 1-D primitive
 *  @param L                the (directional coordinate of the) center of the right 1-D primitive
 *  @param j                the Cartesian exponent of the right 1-D primitive
 * 
 *  @return the linear momentum integral over the two given 1-D primitives
 */
PrimitiveLinearMomentumIntegralEngine::IntegralScalar PrimitiveLinearMomentumIntegralEngine::calculate1D(const double alpha, const double K, const int i, const double beta, const double L, const int j) {

    PrimitiveOverlapIntegralEngine overlap_engine;

    using namespace GQCP::literals;
    return 2.0 * 1.0_ii * beta * overlap_engine.calculate1D(alpha, K, i, beta, L, j + 1) -
           1.0_ii * static_cast<double>(j) * overlap_engine.calculate1D(alpha, K, i, beta, L, j - 1);
}


}  // namespace GQCP
