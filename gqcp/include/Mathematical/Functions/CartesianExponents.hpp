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

#include <array>
#include <cstddef.h>
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
     *  @param array            the array containing the x-, y- and z-exponent (in that order)
     */
    CartesianExponents(const std::array<size_t, 3>& array);

    /**
     *  @param x        the exponent in x
     *  @param y        the exponent in y
     *  @param z        the exponent in z
     */
    CartesianExponents(const size_t x, const size_t y, const size_t z);

    /**
     *  @param vector           the vector containing the x-, y- and z-exponent (in that order) 
     */
    CartesianExponents(const std::vector<size_t>& vector);


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
     *  @return the angular momentum corresponding to these exponents
     */
    size_t angularMomentum() const;

    /**
     *  @return a sorted list of all permutations (i.e. switching x, y, z) of these Cartesian exponents
     */
    std::vector<CartesianExponents> allPermutations() const;

    /**
     *  @return the exponents as an array
     */
    const std::array<size_t, 3>& asArray() const { return this->exponents; }

    /**
     *  @return a textual description of self
     */
    std::string description() const;

    /**
     *  @param direction        the direction (x,y,z) whose exponent should be returned
     *
     *  @return the exponent in the given direction
     */
    size_t value(const CartesianDirection direction) const { return this->exponents[direction]; }
};


}  // namespace GQCP
