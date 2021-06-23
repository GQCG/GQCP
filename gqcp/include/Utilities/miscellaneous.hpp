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


#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iterator>
#include <string>
#include <vector>


namespace GQCP {


/**
 *  A generalization of std::find
 * 
 *  @tparam T           the type of elements stored in the vector
 * 
 *  @param vector       the vector containing the elements
 *  @apram value        the value that should be found
 * 
 *  @return the index of the element that should be found in the given vector
 */
template <typename T>
size_t findElementIndex(const std::vector<T>& vector, const T& value) {

    const auto& it = std::find(vector.begin(), vector.end(), value);  // 'it' for iterator

    // Check if the value was found
    if (it == vector.end()) {
        throw std::out_of_range("findElementIndex(const std::vector<T>&, const T&): the given value was not found in the given vector");
    }

    return std::distance(vector.begin(), it);  // the 'difference' between two iterators is an index
}


/**
 *  Partition a positive integer into its partitions.
 *
 *  @param n        the integer whose partitions are sought
 *  @param k        the partition size
 *
 *  @return a vector of the k-sized partitions
 * 
 *  @example An examples of a 3-way partitions of '2' is:
 *              {2, 0, 0} - {1, 1, 0} - {1, 0, 1} - {0, 2, 0} - {0, 1, 1} - {0, 0, 2}
 */
std::vector<std::vector<size_t>> generatePartitionsOf(const size_t n, const size_t k);

/**
 *  Partition a positive integer into its unique partitions.
 *
 *  @param n        the integer whose partitions are sought
 *  @param k        the partition size
 *
 *  @return the vector of the k-sized partitions
 * 
 *  @example Some examples of unique 3-way partitions a:
 *              {3, 0, 0} - {2, 1, 0} - {1, 1, 1}
 *              {4, 0, 0} - {3, 1, 0} - {2, 2, 0} - {2, 1, 1}
 */
std::vector<std::vector<size_t>> generateUniquePartitionsOf(const size_t n, const size_t k);

/**
 *  @param S    the positive integer to be converted to Gray code
 *
 *  @return the Gray code of the given integer number as a bitset
 */
size_t grayCodeOf(const size_t S);

/**
 *  @param v            the vector index
 *  @param cols         the number of columns in the matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the row-major major (non-contiguous) index given the corresponding vector index
 */
size_t matrixIndexMajor(const size_t v, const size_t cols, const size_t skipped = 0);

/**
 *  @param v            the vector index
 *  @param cols         the number of columns in the matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the row-major minor (contiguous) index given the corresponding vector index
 */
size_t matrixIndexMinor(const size_t v, const size_t cols, const size_t skipped = 0);

/**
 *  Print the time a function takes to be executed
 *
 *  @param method_name      the name of function that is to be executed
 *  @param function         the function call to be made
 */
void printExecutionTime(const std::string& method_name, const std::function<void()>& function);

/**
 *  @param x        the number
 *
 *  @return the strict triangular root of the given number. This is also the dimension of the square matrix whose strict lower/upper triangle has the given number of elements
 */
size_t strictTriangularRootOf(const size_t x);

/**
 *  @param x        the number
 *
 *  @return the triangular root of the given number. This is also the dimension of the square matrix whose lower/upper triangle has the given number of elements
 */
size_t triangularRootOf(const size_t x);

/**
 *  @param filename         the name of the file that should be opened
 *  @param extension        the expected extension of the filename
 */
std::ifstream validateAndOpen(const std::string& filename, const std::string& extension);

/**
 *  @param i            the row index
 *  @param j            the column index
 *  @param cols         the number of columns in de matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the vector index given the corresponding row-major matrix indices
 */
size_t vectorIndex(const size_t i, const size_t j, const size_t cols, size_t skipped = 0);


}  // namespace GQCP
