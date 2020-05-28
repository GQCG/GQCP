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

#include "Utilities/miscellaneous.hpp"

#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>

#include <chrono>
#include <iostream>


namespace GQCP {


/**
 *  @param S    the positive integer to be converted to Gray code
 *
 *  @return the Gray code of the given integer number as a bitset
 */
size_t grayCodeOf(const size_t S) {

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
size_t matrixIndexMajor(const size_t v, const size_t cols, const size_t skipped) {
    return v / (cols - skipped);
}


/**
 *  @param v            the vector index
 *  @param cols         the number of columns in the matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the row-major minor (contiguous) index given the corresponding vector index
 */
size_t matrixIndexMinor(const size_t v, const size_t cols, const size_t skipped) {
    return v % (cols - skipped) + skipped;
}


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
 *  @param x        the number
 *
 *  @return the strict triangular root of the given number. This is also the dimension of the square matrix whose strict lower/upper triangle has the given number of elements
 */
size_t strictTriangularRootOf(const size_t x) {

    return triangularRootOf(x) + 1;
}


/**
 *  @param x        the number
 *
 *  @return the triangular root of the given number. This is also the dimension of the square matrix whose lower/upper triangle has the given number of elements
 */
size_t triangularRootOf(const size_t x) {

    size_t n = static_cast<size_t>((std::sqrt(8 * x + 1) - 1) / 2);

    if (n * (n + 1) != 2 * x) {
        throw std::invalid_argument("triangularRootOf(const size_t): The given number does not have a triangular root.");
    }

    return n;
}


/**
 *  @param i            the row index
 *  @param j            the column index
 *  @param cols         the number of columns in de matrix
 *  @param skipped      the number of columns that are skipped in the matrix representation
 *
 *  @return the vector index given the corresponding row-major matrix indices
 */
size_t vectorIndex(const size_t i, const size_t j, const size_t cols, const size_t skipped) {
    return (j - skipped) + (cols - skipped) * i;
}


/**
 *  @param filename         the name of the file that should be opened
 *  @param extension        the expected extension of the filename
 */
std::ifstream validateAndOpen(const std::string& filename, const std::string& extension) {

    // Find the extension of the given path (https://stackoverflow.com/a/51992)
    std::string filename_extension;  // the extension of the given filename
    std::string::size_type idx = filename.rfind('.');

    if (idx != std::string::npos) {
        filename_extension = filename.substr(idx + 1);
    } else {
        throw std::invalid_argument("validateAndOpen(const std::string&, const std::string&): I did not find an extension in your given file name.");
    }

    if (!(filename_extension == extension)) {
        throw std::invalid_argument("validateAndOpen(const std::string&, const std::string&): The given filen name does not have the expected extension");
    }

    // If the filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
    std::ifstream input_file_stream {filename};
    if (!input_file_stream.good()) {
        throw std::invalid_argument("validateAndOpen(const std::string&, const std::string&): The provided file name is illegible. Maybe you specified a wrong path?");
    } else {
        return input_file_stream;
    }
}


}  // namespace GQCP
