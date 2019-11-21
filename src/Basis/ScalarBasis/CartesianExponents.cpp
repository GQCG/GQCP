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
#include "Basis/ScalarBasis/CartesianExponents.hpp"

#include <algorithm>


namespace GQCP {

/*
 *  CONSTRUCTORS
 */

/**
 *  @param x        the exponent in x
 *  @param y        the exponent in y
 *  @param z        the exponent in z
 */
CartesianExponents::CartesianExponents(size_t x, size_t y, size_t z) :
    x (x),
    y (y),
    z (z)
{}


/**
 *  @param arr      the array containing the x-, y- and z-exponent in that order
 */
CartesianExponents::CartesianExponents(const std::array<size_t, 3>& arr) :
    x (arr[0]),
    y (arr[1]),
    z (arr[2])
{}



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
    if (this->angularMomentum() < rhs.angularMomentum()) { return true; }
    else if (this->angularMomentum() > rhs.angularMomentum()) { return false; }

    else {
        // Compare x-exponents
        if (this->x > rhs.x) { return true; }
        else if (this->x < rhs.x) { return false; }

        else {
            // Compare y-exponents
            if (this->y > rhs.y) { return true; }
            else if (this->y < rhs.y) { return false; }

            else {
                // Compare z-exponents
                if (this->z >= rhs.z) { return true; }  // if, at the end, the z exponents are equal, then the Cartesian exponents should be considered equal
                else { return false; }
            }  // same y
        }  // same x
    }  // same angular momentum
}


/**
 *  @param rhs      the right-hand side of the operator ==
 *
 *  @return if the Cartesian exponents are considered equal
 */
bool CartesianExponents::operator==(const CartesianExponents& rhs) const {
    return (this->x == rhs.x) && (this->y == rhs.y) && (this->z == rhs.z);
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
 *  @param direction        the direction (x,y,z) whose exponent should be returned
 *
 *  @return the exponent in the given direction
 */
size_t CartesianExponents::value(CartesianDirection direction) const {

    switch (direction) {
        case CartesianDirection::x:
            return this->x;

        case CartesianDirection::y:
            return this->y;

        case CartesianDirection::z:
            return this->z;
    }
}


/**
 *  @return the angular momentum corresponding to these exponents
 */
size_t CartesianExponents::angularMomentum() const {
    return this->x + this->y + this->z;
}


/**
 *  @return the exponents as an array
 */
std::array<size_t, 3> CartesianExponents::asArray() const {
    return std::array<size_t, 3> {this->x, this->y, this->z};
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
