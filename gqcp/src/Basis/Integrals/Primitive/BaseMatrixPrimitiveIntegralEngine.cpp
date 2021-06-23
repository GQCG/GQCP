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

#include "Basis/Integrals/Primitive/BaseMatrixPrimitiveIntegralEngine.hpp"

#include <stdexcept>


namespace GQCP {


/*
 *  MARK: Constructors
 */

/** 
 *  @param component                    The initial component of the dyadic Cartesian operator over which the primitive engine is prepared to calculate integrals.
 * 
 *  @note The primitive engine's state may be changed through the `prepareStateForComponent` method.
 */
BaseMatrixPrimitiveIntegralEngine::BaseMatrixPrimitiveIntegralEngine(const DyadicCartesianDirection component) :
    component {component} {}


/*
 *  MARK: Destructor
 */

/**
 *  A pure virtual destructor in order to make this class abstract.
 */
BaseMatrixPrimitiveIntegralEngine::~BaseMatrixPrimitiveIntegralEngine() {}


/*
 *  MARK: Components
 */

/**
 *  Prepare this engine's internal state such that it is able to calculate integrals over the given component of the operator.
 * 
 *  @param component                The index of the component of the operator.
 * 
 *  @note See also `DyadicCartesianDirection`.
 */
void BaseMatrixPrimitiveIntegralEngine::prepareStateForComponent(const size_t component) {

    switch (component) {
    case 0: {
        this->component = DyadicCartesianDirection::xx;
        break;
    }

    case 1: {
        this->component = DyadicCartesianDirection::xy;
        break;
    }

    case 2: {
        this->component = DyadicCartesianDirection::xz;
        break;
    }

    case 3: {
        this->component = DyadicCartesianDirection::yx;
        break;
    }

    case 4: {
        this->component = DyadicCartesianDirection::yy;
        break;
    }

    case 5: {
        this->component = DyadicCartesianDirection::yz;
        break;
    }

    case 6: {
        this->component = DyadicCartesianDirection::zx;
        break;
    }

    case 7: {
        this->component = DyadicCartesianDirection::zy;
        break;
    }

    case 8: {
        this->component = DyadicCartesianDirection::zz;
        break;
    }

    default: {
        throw std::invalid_argument("BaseMatrixPrimitiveIntegralEngine::prepareStateForComponent(const size_t): The given component is out of bound for the number of components of a matrix primitive integral engine.");
        break;
    }
    }
}


}  // namespace GQCP
