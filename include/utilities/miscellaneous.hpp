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
#ifndef GQCP_MISCELLANEOUS_HPP
#define GQCP_MISCELLANEOUS_HPP


#include <string>
#include <functional>

#include <stdlib.h>



#include <iostream>



namespace GQCP {


/**
 *  Print the time a function takes to be executed
 *
 *  @param method_name      the name of function that is to be executed
 *  @param function         the function call to be made
 */
void printExecutionTime(const std::string& method_name, const std::function<void()>& function);


/**
 *  @param S    the positive integer to be converted to Gray code
 *
 *  @return the Gray code of the given integer number as a bitset
 */
size_t gray_code(size_t S);

/**
 *  @param v            the vector index
 *  @param cols         the number of columns in the matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the row-major major (non-contiguous) index given the corresponding vector index
 */
size_t matrixIndexMajor(size_t v, size_t cols, size_t skipped=0);

/**
 *  @param v            the vector index
 *  @param cols         the number of columns in the matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the row-major minor (contiguous) index given the corresponding vector index
 */
size_t matrixIndexMinor(size_t v, size_t cols, size_t skipped=0);

/**
 *  @param i            the row index
 *  @param j            the column index
 *  @param cols         the number of columns in de matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the vector index given the corresponding row-major matrix indices
 */
size_t vectorIndex(size_t i, size_t j, size_t cols, size_t skipped=0);
    

/**
 *  Partition a positive integer into its unique partitions
 *
 *  @tparam k       the partition size
 *
 *  @param n        the integer whose partitions are sought
 *
 *  @return the vector of the k-sized partitions
 */
template <size_t k>
std::vector<std::array<size_t, k>> uniquePartitions(size_t n) {

    static_assert(k > 0, "template<size_t> uniquePartitions(size_t): the template parameter must be larger than zero");


    // The main algorithm starts from {n, 0, ..., 0} and moves a 1 from the right-most number (>1) to the left-most position that holds a value at least 2 smaller. If there are none such numbers left, the algorithm is finished
    // Note that:
    //  - The largest number L in the partition will always be at the first position; the partition will always be sorted from largest to smallest value
    //  - The right-most 1s can be ignored, since moving a 1 cannot create a new partition
    //  - Moving a 1 to a position that holds a value L-1 or higher can be ignored, since it will always create a permutation of a previous partition
    // Some examples of 3-way partitions are:
    //  {3, 0, 0} - {2, 1, 0} - {1, 1, 1}
    //  {4, 0, 0} - {3, 1, 0} - {2, 2, 0} - {2, 1, 1}


    std::array<size_t, k> partition {};
    partition[0] = n;  // start with {n, 0, ..., 0}
    std::vector<std::array<size_t, k>> unique_partitions {partition};
    while (true) {

        // Find the right-most position that holds a value larger than 1
        auto subtraction_it = std::find_if(partition.rbegin(), partition.rend(), [](size_t x) { return x > 1; });
        if (subtraction_it == partition.rend()) {
            break;
        }

        // Find the left-most position that is smaller than (value - 1)
        size_t value = *subtraction_it;
        auto addition_it = std::find_if(partition.begin(), partition.end(), [value](size_t x) { return x < value - 1; });
        if (addition_it == partition.end()) {
            break;
        }

        // If there are such positions, proceed to move a 1
        (*subtraction_it)--;
        (*addition_it)++;
        unique_partitions.push_back(partition);
    };

    return unique_partitions;
}


}  // namespace GQCP


#endif  // GQCP_MISCELLANEOUS_HPP
