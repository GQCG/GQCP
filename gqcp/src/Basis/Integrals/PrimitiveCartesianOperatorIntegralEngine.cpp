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

#include "Basis/Integrals/PrimitiveCartesianOperatorIntegralEngine.hpp"

#include <stdexcept>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Construct a PrimitiveCartesianOperatorIntegralEngine from its members.
 * 
 *  @param component                    the initial component of the Cartesian operator this engine should calculate integrals over
 * 
 *  @note The primitive engine's state may be changed through the prepareStateForComponent method.
 */
PrimitiveCartesianOperatorIntegralEngine::PrimitiveCartesianOperatorIntegralEngine(const CartesianDirection component) :
    component {component} {}


/*
 *  DESTRUCTOR
 */

/**
 *  A pure virtual destructor should have an empty body.
 */
PrimitiveCartesianOperatorIntegralEngine::~PrimitiveCartesianOperatorIntegralEngine() {}


/*
 *  PUBLIC METHODS
 */

/**
 *  Prepare this engine's internal state such that it is able to calculate integrals over the given component of the operator.
 * 
 *  @param component                the index of the component of the operator
 */
void PrimitiveCartesianOperatorIntegralEngine::prepareStateForComponent(const size_t component) {

    switch (component) {
    case 0: {
        this->component = CartesianDirection::x;
        break;
    }

    case 1: {
        this->component = CartesianDirection::y;
        break;
    }

    case 2: {
        this->component = CartesianDirection::z;
        break;
    }

    default: {
        throw std::invalid_argument("PrimitiveDipoleIntegralEngine::prepareStateForComponent(const size_t): The given component is out of bound for the number of components of the primitive dipole integral engine.");
        break;
    }
    }
}


}  // namespace GQCP
