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

#include "Basis/Integrals/PrimitiveDipoleIntegralEngine.hpp"

#include "Basis/Integrals/McMurchieDavidsonCoefficient.hpp"
#include "Basis/Integrals/PrimitiveOverlapIntegralEngine.hpp"

#include <boost/math/constants/constants.hpp>


namespace GQCP {


/**
 *  Construct a PrimitiveDipoleIntegralEngine from its members.
 * 
 *  @param dipole_operator              the dipole operator over which this engine should calculate integrals
 *  @param component                    the current component of the dipole operator this engine can calculate integrals over
 */
PrimitiveDipoleIntegralEngine::PrimitiveDipoleIntegralEngine(const ElectronicDipoleOperator& dipole_operator, const CartesianDirection component) :
    dipole_operator {dipole_operator},
    PrimitiveCartesianOperatorIntegralEngine(component) {}


// PUBLIC METHODS

/**
 *  @param left             the left Cartesian GTO (primitive)
 *  @param right            the right Cartesian GTO (primitive)
 * 
 *  @return the dipole integral (of the current component) over the two given primitives
 */
PrimitiveDipoleIntegralEngine::IntegralScalar PrimitiveDipoleIntegralEngine::calculate(const CartesianGTO& left, const CartesianGTO& right) {

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

    PrimitiveOverlapIntegralEngine<GTOShell> overlap_engine;


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
 *  @return the dipole integral over the two given 1-D primitives
 */
PrimitiveDipoleIntegralEngine::IntegralScalar PrimitiveDipoleIntegralEngine::calculate1D(const double alpha, const double K, const int i, const double beta, const double L, const int j) {

    // Prepare some variables.
    const auto p = alpha + beta;
    const auto P = (alpha * K + beta * L) / p;  // one of the components of the center of mass of the Gaussian overlap distribution

    const auto Delta_PO = P - this->dipole_operator.reference()(this->component);  // one of the components of the distance of P and the origin of the dipole operator

    // Calculate the dipole integral over the current component.
    const McMurchieDavidsonCoefficient E {K, alpha, L, beta};
    return -std::pow(boost::math::constants::pi<IntegralScalar>() / p, 0.5) * (E(i, j, 1) + Delta_PO * E(i, j, 0));  // the minus sign comes from the electronic dipole operator
}


}  // namespace GQCP
