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
    std::sort(index_vector.begin(), index_vector.end());
    auto iterator = std::unique(index_vector.begin(), index_vector.end());

    // If nothing was removed by `std::unique`, i.e. all indices were unique, the return value of the iterator shall equal the last value of the vector.
    // Since the operator string will not be zero in this case, we must return false.
    return (iterator != index_vector.end());
}


/*
*  MARK: Public methods
*/

int SpinUnresolvedOperatorString::phaseFactorAfterSorting() {

    auto index_vector = this->operatorIndices();
    const auto number_of_occupied_spinors = index_vector.size();
    int phase_factor = 1;

    for (int i = 0; i < number_of_occupied_spinors; ++i) {

        for (int j = i; j >= 0; --j) {
            if (index_vector[j] > index_vector[i]) {
                phase_factor *= -1;
            }
        }
    }
    return phase_factor;
}


void SpinUnresolvedOperatorString::sort() {

    const auto phase_factor = this->phaseFactorAfterSorting();

    std::sort(this->indices.begin(), this->indices.end());
    this->p *= phase_factor;
}


std::vector<SpinUnresolvedOperatorString> SpinUnresolvedOperatorString::splitIntoSystemAndEnvironment(const std::vector<char>& partition) {

    const auto& index_vector = this->operatorIndices();

    std::vector<size_t> index_vector_I;
    std::vector<size_t> index_vector_J;
    int phase_factor = this->phaseFactor();

    for (int i = 0; i < index_vector.size(); ++i) {

        if (partition[i] == 'I') {
            for (int j = i; j >= 0; --j) {
                if (partition[j] == 'J') {
                    phase_factor *= -1;
                }
            }
            index_vector_I.push_back(index_vector[i]);
        } else {
            index_vector_J.push_back(index_vector[i]);
        }
    }

    const SpinUnresolvedOperatorString operator_string_I {index_vector_I, phase_factor};
    const SpinUnresolvedOperatorString operator_string_J {index_vector_J, 1};
    return {operator_string_I, operator_string_J};
}


}  // namespace GQCP