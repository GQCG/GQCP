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

#include "Basis/Integrals/PrimitiveOverlapIntegralEngine.hpp"

#include "Basis/Integrals/McMurchieDavidsonCoefficient.hpp"

#include <boost/math/constants/constants.hpp>


namespace GQCP {


/*
 *  PUBLIC METHODS
 */

/**
 *  @param left             the left Cartesian GTO (primitive)
 *  @param right            the right Cartesian GTO (primitive)
 * 
 *  @return the overlap integral over the two given primitives
 */
double PrimitiveOverlapIntegralEngine::calculate(const CartesianGTO& left, const CartesianGTO& right) {

    double primitive_integral = 1.0;
    for (const auto& direction : {GQCP::CartesianDirection::x, GQCP::CartesianDirection::y, GQCP::CartesianDirection::z}) {
        const auto i = left.cartesianExponents().value(direction);
        const auto j = right.cartesianExponents().value(direction);

        primitive_integral *= this->calculate1D(left.gaussianExponent(), left.center()(direction), i, right.gaussianExponent(), right.center()(direction), j);
    }

    return primitive_integral;
}


/**
 *  @param alpha            the Gaussian exponent of the left 1-D primitive
 *  @param K                the (directional coordinate of the) center of the left 1-D primitive
 *  @param i                the Cartesian exponent of the left 1-D primitive
 *  @param beta             the Gaussian exponent of the right 1-D primitive
 *  @param L                the (directional coordinate of the) center of the right 1-D primitive
 *  @param j                the Cartesian exponent of the right 1-D primitive
 * 
 *  @return the overlap integral over the two given 1-D primitives
 */
double PrimitiveOverlapIntegralEngine::calculate1D(const double alpha, const double K, const int i, const double beta, const double L, const int j) {

    const auto p = alpha + beta;
    const McMurchieDavidsonCoefficient E {K, alpha, L, beta};

    return std::pow(boost::math::constants::pi<double>() / p, 0.5) * E(i, j, 0);
}


}  // namespace GQCP
