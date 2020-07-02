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

#include "Basis/Integrals/PrimitiveAngularMomentumIntegralEngine.hpp"

#include "Basis/Integrals/PrimitiveDipoleIntegralEngine.hpp"
#include "Basis/Integrals/PrimitiveLinearMomentumIntegralEngine.hpp"
#include "Basis/Integrals/PrimitiveOverlapIntegralEngine.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Construct a PrimitiveAngularMomentumIntegralEngine from its members.
 * 
 *  @param angular_momentum_operator                the angular momentum operator over which the integrals should be calculated
 *  @param component                                the initial component of the angular momentum operator this engine should calculate integrals over
 */
PrimitiveAngularMomentumIntegralEngine::PrimitiveAngularMomentumIntegralEngine(const AngularMomentumOperator& angular_momentum_operator, const CartesianDirection component) :
    angular_momentum_operator {angular_momentum_operator},
    PrimitiveCartesianOperatorIntegralEngine(component) {}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param left             the left Cartesian GTO (primitive)
 *  @param right            the right Cartesian GTO (primitive)
 * 
 *  @return the angular momentum integral (of the current component) over the two given primitives
 */
PrimitiveAngularMomentumIntegralEngine::IntegralScalar PrimitiveAngularMomentumIntegralEngine::calculate(const CartesianGTO& left, const CartesianGTO& right) {

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


    // For each component of the angular momentum operator, the integrals can be calculated through overlap integrals, linear momentum integrals and position/dipole integrals.
    PrimitiveOverlapIntegralEngine overlap_engine;
    PrimitiveLinearMomentumIntegralEngine linear_momentum_engine;
    PrimitiveDipoleIntegralEngine dipole_engine {ElectronicDipoleOperator(this->angular_momentum_operator.reference())};


    // We'll have to switch sign because we use dipole integrals instead of position integrals.
    switch (this->component) {
    case CartesianDirection::x: {
        dipole_engine.prepareStateForComponent(CartesianDirection::y);
        linear_momentum_engine.prepareStateForComponent(CartesianDirection::z);
        const IntegralScalar term1 = -dipole_engine.calculate1D(alpha, K_y, k, beta, L_y, l) * linear_momentum_engine.calculate1D(alpha, K_z, m, beta, L_z, n);  // incorporate the sign

        linear_momentum_engine.prepareStateForComponent(CartesianDirection::y);
        dipole_engine.prepareStateForComponent(CartesianDirection::z);
        const IntegralScalar term2 = -linear_momentum_engine.calculate1D(alpha, K_y, k, beta, L_y, l) * dipole_engine.calculate1D(alpha, K_z, m, beta, L_z, n);  // incorporate the sign

        return overlap_engine.calculate1D(alpha, K_x, i, beta, L_x, j) * (term1 - term2);  // the cross product
        break;
    }

    case CartesianDirection::y: {
        dipole_engine.prepareStateForComponent(CartesianDirection::z);
        linear_momentum_engine.prepareStateForComponent(CartesianDirection::x);
        const IntegralScalar term1 = -dipole_engine.calculate1D(alpha, K_z, m, beta, L_z, n) * linear_momentum_engine.calculate1D(alpha, K_x, i, beta, L_x, j);  // incorporate the sign

        linear_momentum_engine.prepareStateForComponent(CartesianDirection::z);
        dipole_engine.prepareStateForComponent(CartesianDirection::x);
        const IntegralScalar term2 = -linear_momentum_engine.calculate1D(alpha, K_z, m, beta, L_z, n) * dipole_engine.calculate1D(alpha, K_x, i, beta, L_x, j);  // incorporate the sign

        return overlap_engine.calculate1D(alpha, K_y, k, beta, L_y, l) * (term1 - term2);  // the cross product
        break;
    }

    case CartesianDirection::z: {
        dipole_engine.prepareStateForComponent(CartesianDirection::x);
        linear_momentum_engine.prepareStateForComponent(CartesianDirection::y);
        const IntegralScalar term1 = -dipole_engine.calculate1D(alpha, K_x, i, beta, L_x, j) * linear_momentum_engine.calculate1D(alpha, K_y, k, beta, L_y, l);  // incorporate the sign

        linear_momentum_engine.prepareStateForComponent(CartesianDirection::x);
        dipole_engine.prepareStateForComponent(CartesianDirection::y);
        const IntegralScalar term2 = -linear_momentum_engine.calculate1D(alpha, K_x, i, beta, L_x, j) * dipole_engine.calculate1D(alpha, K_y, k, beta, L_y, l);  // incorporate the sign

        return overlap_engine.calculate1D(alpha, K_z, m, beta, L_z, n) * (term1 - term2);  // the cross product
        break;
    }
    }
}


}  // namespace GQCP
