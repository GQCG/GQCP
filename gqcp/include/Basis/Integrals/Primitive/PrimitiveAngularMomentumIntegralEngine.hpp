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

#include "Basis/Integrals/Primitive/BaseVectorPrimitiveIntegralEngine.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Operator/FirstQuantized/AngularMomentumOperator.hpp"


namespace GQCP {


/**
 *  A one-electron primitive integral engine that can calculate integrals over the angular momentum operator.
 */
class PrimitiveAngularMomentumIntegralEngine:
    public BaseVectorPrimitiveIntegralEngine {
public:
    // The scalar representation of an integral.
    using IntegralScalar = AngularMomentumOperator::Scalar;

    // The type of shell that this integral engine is related to.
    using Shell = GTOShell;

    // The type of primitive that underlies the type of shell.
    using Primitive = Shell::Primitive;


private:
    AngularMomentumOperator angular_momentum_operator;  // the angular momentum operator over which the integrals should be calculated


public:
    // CONSTRUCTORS

    /**
     *  Construct a PrimitiveAngularMomentumIntegralEngine from its members.
     * 
     *  @param angular_momentum_operator                the angular momentum operator over which the integrals should be calculated
     *  @param component                                the initial component of the angular momentum operator this engine should calculate integrals over
     */
    PrimitiveAngularMomentumIntegralEngine(const AngularMomentumOperator& angular_momentum_operator, const CartesianDirection component = CartesianDirection::x);


    // PUBLIC METHODS

    /**
     *  @param left             the left Cartesian GTO (primitive)
     *  @param right            the right Cartesian GTO (primitive)
     * 
     *  @return the angular momentum integral (of the current component) over the two given primitives
     */
    IntegralScalar calculate(const CartesianGTO& left, const CartesianGTO& right);
};


}  // namespace GQCP
