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

#include "ONVBasis/SpinUnresolvedOperatorString.hpp"


namespace GQCP {


/*
 *  MARK: General information
 */

/**
 *  Check whether the operator string in question will result in zero when applied to the wave function.
 *
 *  Note: There are different cases when an operator string will result in a zero value. This method checks all of them.
 */
bool SpinUnresolvedOperatorString::isZero() const {

    // Check if the same index appears more than once in the operator string.
    // Start by initializing a copy of the index vector associated with the operator string.
    auto index_vector = this->operatorIndices();

    // Since the operator string represents either only annihilation or creation operators, repetition of an index means that that index will be annihilated or created twice in any ONV following the operator string, which will automatically result in zero.
    // We will use the std::unique function to check this condition. `std::unique` needs a sorted vector and removes all but the first instance of any unique group of elements. It then returns a sequence that's not necessarily equal to the original vector, as duplicate values are removed. Full explanation can be found here: https://stackoverflow.com/questions/46477764/check-stdvector-has-duplicates/46477901.
    sort(index_vector.begin(), index_vector.end());
    auto iterator = std::unique(index_vector.begin(), index_vector.end());

    // If nothing was removed by `std::unique`, i.e. all indices were unique, the return value of the iterator shall equal the last value of the vector.
    // Since the operator string will not be zero in this case, we must return false.
    return (iterator != index_vector.end());
}


}  // namespace GQCP