// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2019  the GQCG developers
//
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
//
#pragma once


#include "Mathematical/CartesianDirection.hpp"

#include <array>
#include <vector>


namespace GQCP {


/**
 *  A class that represents exponents of the Cartesian functions x, y and z
 */
struct CartesianExponents {
public:
    std::array<size_t, 3> exponents;  // array containing the x, y, z exponents (in that order)


public:
    // CONSTRUCTORS
    /**
     *  @param x        the exponent in x
     *  @param y        the exponent in y
     *  @param z        the exponent in z
     */
    CartesianExponents(size_t x, size_t y, size_t z);

    /**
     *  @param arr      the array containing the x-, y- and z-exponent in that order
     */
    CartesianExponents(const std::array<size_t, 3>& arr);


    // OPERATORS
    /**
     *  @param rhs      the right-hand side of the operator <
     *
     *  @return if these Cartesian exponents are 'smaller' than the ones on the right-hand side. The following logic is used: lhs < rhs
     *      - if lhs's angular momentum is smaller
     *      - if both angular momenta are equal, x takes precedence over y, over z
     *
     *  This means that {1, 0, 0}(=x) < {2, 0, 0}(=x^2), and {2, 0, 0}(=x^2) < {1, 1, 0}(=xy)
     */
    bool operator<(const CartesianExponents& rhs) const;

    /**
     *  @param rhs      the right-hand side of the operator ==
     *
     *  @return if the Cartesian exponents are considered equal
     */
    bool operator==(const CartesianExponents& rhs) const;

    /**
     *  @param rhs      the right-hand side of the operator ==
     *
     *  @return if the Cartesian exponents are considered different
     */
    bool operator!=(const CartesianExponents& rhs) const;


    // PUBLIC METHODS
    /**
     *  @param direction        the direction (x,y,z) whose exponent should be returned
     *
     *  @return the exponent in the given direction
     */
    size_t value(CartesianDirection direction) const;

    /**
     *  @return the angular momentum corresponding to these exponents
     */
    size_t angularMomentum() const;

    /**
     *  @return the exponents as an array
     */
    const std::array<size_t, 3>& asArray() const;

    /**
     *  @return a sorted list of all permutations (i.e. switching x, y, z) of these Cartesian exponents
     */
    std::vector<CartesianExponents> allPermutations() const;
};


}  // namespace GQCP
