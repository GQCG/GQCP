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


#include "Mathematical/Functions/CartesianDirection.hpp"

#include <cstddef>


namespace GQCP {


/**
 *  A base class for integral engines over an operator with three vector components (x, y, z).
 */
class BaseVectorPrimitiveIntegralEngine {
public:
protected:
    // The component over which the primitive engine is prepared to calculate integrals.
    CartesianDirection component;


public:
    /*
     *  MARK: Constructors
     */

    /** 
     *  @param component                    The initial component of the Cartesian operator over which the primitive engine is prepared to calculate integrals.
     * 
     *  @note The primitive engine's state may be changed through the `prepareStateForComponent` method.
     */
    BaseVectorPrimitiveIntegralEngine(const CartesianDirection component = CartesianDirection::x);


    /*
     *  MARK: Destructor
     */

    /**
     *  A pure virtual destructor in order to make this class abstract.
     */
    virtual ~BaseVectorPrimitiveIntegralEngine() = 0;


    /*
     *  MARK: Components
     */

    // The number of components the operator has. For a vector operator, this is equal to 3.
    static constexpr size_t Components = 3;

    /**
     *  Prepare this engine's internal state such that it is able to calculate integrals over the given component of the operator.
     * 
     *  @param component                The index of the component of the operator.
     * 
     *  @note See also `CartesianDirection`.
     */
    void prepareStateForComponent(const size_t component);
};


}  // namespace GQCP
