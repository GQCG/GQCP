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

#include "Basis/ScalarBasis/CartesianExponents.hpp"

#include <algorithm>


namespace GQCP {

/*
 *  CONSTRUCTORS
 */

/**
 *  @param arr      the array containing the x-, y- and z-exponent in that order
 */
CartesianExponents::CartesianExponents(const std::array<size_t, 3>& arr) :
    exponents {arr} {}


/**
 *  @param x        the exponent in x
 *  @param y        the exponent in y
 *  @param z        the exponent in z
 */
CartesianExponents::CartesianExponents(const size_t x, const size_t y, const size_t z) :
    exponents {x, y, z} {}


/*
 *  OPERATORS
 */

/**
 *  @param rhs      the right-hand side of the operator <
 *
 *  @return if these Cartesian exponents are 'smaller' than the ones on the right-hand side. The following logic is used: lhs < rhs
 *      - if lhs's angular momentum is smaller
 *      - if both angular momenta are equal, x takes precedence over y, over z
 *
 *  This means that {1, 0, 0}(=x) < {2, 0, 0}(=x^2), and {2, 0, 0}(=x^2) < {1, 1, 0}(=xy)
 */
bool CartesianExponents::operator<(const CartesianExponents& rhs) const {

    // Compare angular momenta
    if (this->angularMomentum() < rhs.angularMomentum()) {
        return true;
    } else if (this->angularMomentum() > rhs.angularMomentum()) {
        return false;
    } else {
        // Compare all exponents for x -> y -> z, x takes precedence over y, over z
        for (const auto& direction : {CartesianDirection::x, CartesianDirection::y, CartesianDirection::z}) {
            if (this->exponents[direction] > rhs.value(direction)) {
                return true;
            }
            if (this->exponents[direction] < rhs.value(direction)) {
                return false;
            }
        }  // if, at the end, all exponents are equal, then the Cartesian exponents should be considered equal
        return false;
    }  // same angular momentum
}


/**
 *  @param rhs      the right-hand side of the operator ==
 *
 *  @return if the Cartesian exponents are considered equal
 */
bool CartesianExponents::operator==(const CartesianExponents& rhs) const {
    return (this->value(CartesianDirection::x) == rhs.value(CartesianDirection::x)) && (this->value(CartesianDirection::y) == rhs.value(CartesianDirection::y)) && (this->value(CartesianDirection::z) == rhs.value(CartesianDirection::z));
}


/**
 *  @param rhs      the right-hand side of the operator ==
 *
 *  @return if the Cartesian exponents are considered different
 */
bool CartesianExponents::operator!=(const CartesianExponents& rhs) const {
    return !(this->operator==(rhs));
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the angular momentum corresponding to these exponents
 */
size_t CartesianExponents::angularMomentum() const {

    size_t momentum = 0;
    for (const auto& direction : {CartesianDirection::x, CartesianDirection::y, CartesianDirection::z}) {
        momentum += this->exponents[direction];
    }
    return momentum;
}


/**
 *  @return a sorted list of all permutations (i.e. switching x, y, z) of these Cartesian exponents
 */
std::vector<CartesianExponents> CartesianExponents::allPermutations() const {

    std::vector<CartesianExponents> all_permutations;
    all_permutations.reserve(6);  // maximally 6 = 3! possible permutations

    // Generate all permutations using this' array representation: sort and then next_permutation
    auto arr = this->asArray();  // arr: array
    std::sort(arr.begin(), arr.end());
    do {
        all_permutations.emplace_back(arr);
    } while (std::next_permutation(arr.begin(), arr.end()));


    std::sort(all_permutations.begin(), all_permutations.end());
    return all_permutations;
}


}  // namespace GQCP
