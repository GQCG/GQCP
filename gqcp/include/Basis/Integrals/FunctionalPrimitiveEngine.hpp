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

#include "Mathematical/Functions/CartesianGTO.hpp"

#include <functional>


namespace GQCP {


/**
 *  A custom primitive engine that can calculate one-electron integrals over Cartesian GTOs according to a custom implementation. It is limited to calculating integrals over scalar operators.
 * 
 *  @param _IntegralScalar          The scalar representation of one of the primitive integrals.
 */
template <typename _IntegralScalar>
class FunctionalPrimitiveEngine {
public:
    // The scalar representation of one of the primitive integrals.
    using IntegralScalar = _IntegralScalar;

    // The number of components of the operator over which the primitive engine can calculate integrals.
    static constexpr auto Components = 1;

private:
    // A user-supplied custom function that can calculate primitive integrals over two Cartesian GTOs.
    std::function<IntegralScalar(const CartesianGTO&, const CartesianGTO&)> function;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param function             A user-supplied custom function that can calculate primitive integrals over two Cartesian GTOs.
     */
    FunctionalPrimitiveEngine(const std::function<IntegralScalar(const CartesianGTO&, const CartesianGTO&)>& function) :
        function {function} {}


    /*
     *  MARK: Integral calculation
     */

    /**
     *  Prepare this engine's internal state such that it is able to calculate integrals over the given component of the operator.
     * 
     *  @param component                The index of the component of the operator.
     * 
     *  @note Since this kind of primitive engine is limited to 1-component operators (i.e. scalar operators), this method has no effect.
     */
    void prepareStateForComponent(const size_t component) {}

    /**
     *  Calculate the primitive integral over two `CartesianGTO`s.
     * 
     *  @param left             The Cartesian GTO of the bra side.
     *  @param right            The Cartesian GTO of the ket site.
     * 
     *  @return The value of the integral over the two given `CartesianGTO`s.
     */
    IntegralScalar calculate(const CartesianGTO& left, const CartesianGTO& right) const { return this->function(left, right); }
};


}  // namespace GQCP
