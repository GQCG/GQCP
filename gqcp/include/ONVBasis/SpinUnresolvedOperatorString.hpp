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


namespace GQCP {


/**
 *  A spin-unresolved operator string.

 *  An spin-unresolved operator string represents a string of either annihilation or creation operators by its indices.
 *  For example, an operator string represented by indices <1, 2, 3> represents either:
 *      a_1^\dagger a_2^\dagger a_3^\dagger
 *  or
 *      a_1 a_2 a_3.
 *  
 *  An operator string is always represented by the indices of the operators. Whether it denotes annihilation or creation operators depends on the context in which an operator string is used.
 *  The operator strings are different from ONV's. Instead of representing the way orbitals are occupied, they purely represent the order of certain operators.
 */
class SpinUnresolvedOperatorString {
private:
    // The vector representing the indices of the annihilation or creation operators in the operator string.
    std::vector<size_t> indices;

public:
    /*
     * MARK: Constructors
     */

    /**
     *  Construct a `SpinUnresolvedOperatorString` from the vector of indices that it encapsulates.
     * 
     *  @param index_vector                The vector containing the operator indices.
     */
    SpinUnresolvedOperatorString(const std::vector<size_t>& index_vector) :
        indices {index_vector} {};


    /*
     *  MARK: General information
     */

    /**
     *  Retrieve the operator indices from the `SpinUnresolvedOPeratorString`.
     */
    std::vector<size_t> operatorIndices() const { return this->indices; }

    /**
     * Check whether the operator string in question will result in zero when applied to the wave function.
     * 
     * Note: There are different cases when an operator string will result in a zero value. This method checks all of them.
     */
    bool isZero() const {

        // Check if adjacent index pairs in the operator string are the same. If there are identical adjacent pairs, the operator stringg will always result in zero due to the fermion anti-commutation rules.
        for (auto left = indices.begin(), right = left + 1, last = indices.end(); right != last; ++left, ++right) {
            if (*left == *right) {
                return true;
            }
        }

        //
    }
};


}  // namespace GQCP
