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
#include "Utilities/miscellaneous.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>

#include <chrono>
#include <iostream>


namespace GQCP {


/**
 *  Print the time a function takes to be executed
 *
 *  @param method_name      the name of function that is to be executed
 *  @param function         the function call to be made
 */
void printExecutionTime(const std::string& method_name, const std::function<void()>& function) {

    // High resolution clock example from (https://stackoverflow.com/a/12231232/7930415)
    auto start = std::chrono::high_resolution_clock::now();

    function();

    auto stop = std::chrono::high_resolution_clock::now();


    // Print the timings
    std::cout << method_name << " took "
              << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
              << " microseconds ("
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
              << " milliseconds) to complete." << std::endl;
}


/**
 *  @param S    the positive integer to be converted to Gray code
 *
 *  @return the Gray code of the given integer number as a bitset
 */
size_t gray_code(size_t S) {

    // See (https://en.wikipedia.org/wiki/Gray_code#Converting_to_and_from_Gray_code)
    return S ^ (S >> 1);
}


/**
 *  @param v            the vector index
 *  @param cols         the number of columns in the matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the row-major major (non-contiguous) index given the corresponding vector index
 */
size_t matrixIndexMajor(size_t v, size_t cols, size_t skipped) {
    return v / (cols - skipped);
}


/**
 *  @param v            the vector index
 *  @param cols         the number of columns in the matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the row-major minor (contiguous) index given the corresponding vector index
 */
size_t matrixIndexMinor(size_t v, size_t cols, size_t skipped) {
    return v % (cols - skipped) + skipped;
}


/**
 *  @param i            the row index
 *  @param j            the column index
 *  @param cols         the number of columns in de matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the vector index given the corresponding row-major matrix indices
 */
size_t vectorIndex(size_t i, size_t j, size_t cols, size_t skipped) {
    return (j - skipped) + (cols - skipped) * i;
}


/**
 *  @param x        the number
 *
 *  @return the triangular root of the given number. This is also the dimension of the square matrix whose lower/upper triangle has the given number of elements
 */
size_t triangularRoot(const size_t x) {

    size_t n = static_cast<size_t>( (std::sqrt(8*x+1) - 1)/2 );

    if (n * (n+1) != 2*x) {
        throw std::invalid_argument("triangularRoot(const size_t): The given number does not have a triangular root.");
    }

    return n;
}


/**
 *  @param x        the number
 *
 *  @return the strict triangular root of the given number. This is also the dimension of the square matrix whose strict lower/upper triangle has the given number of elements
 */
size_t strictTriangularRoot(const size_t x) {

    return triangularRoot(x) + 1;
}


}  // namespace GQCP
