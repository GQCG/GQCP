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
std::vector<std::vector<size_t>> generatePartitionsOf(const size_t n, const size_t k) {

    // Let's take an example generatePartitionsOf(3, 3) to explain why we can use a recursion algorithm. Here are the 3-way partitions of 3:
    //  3 0 0       this is the standard case
    //
    //  2 1 0       this is of the form {2, [2-way partitions of 1]}
    //  2 0 1
    //
    //  1 2 0       this is of the form {1, [2-way partitions of 2]}
    //  1 1 1
    //  1 0 2
    //
    //  0 3 0       this is of the form {0, [2-way partitions of 3]}
    //  0 2 1
    //  0 1 2
    //  0 0 3


    // A base recursion case for k == 1 or k == 0.
    if (k <= 1) {
        return {{n}};
    }

    // A base recursion case for n == 1:  a k-way partition of 1.
    if (n == 1) {
        std::vector<std::vector<size_t>> base_partitions;

        // A k-way partition of 1 is like
        //      {1, 0, 0}
        //      {0, 1, 0}
        //      {0, 0, 1}
        for (size_t i = 0; i < k; i++) {
            std::vector<size_t> partition(k, 0);  // zero-initialize a size-k vector
            partition[i] = 1;

            base_partitions.push_back(partition);
        }

        return base_partitions;
    }


    // If we're here, this means n >= 2 and the actual recursion should start.
    std::vector<size_t> partition(k, 0);

    partition[0] = n;  // start with {n, 0, ..., 0}
    std::vector<std::vector<size_t>> partitions {partition};

    for (int i = n - 1; i >= 0; i--) {

        // Update the leading number.
        partition[0] = i;

        // We must now collect all sub-partitions of size (k-1) of the number (n-i).
        const auto sub_partitions = generatePartitionsOf(n - i, k - 1);

        // And then we append all partitions of the form {i-1, [subpartition]}.
        for (const auto& sub_partition : sub_partitions) {
            std::copy_n(sub_partition.begin(), k - 1, partition.begin() + 1);

            partitions.push_back(partition);
        }
    }

    return partitions;
}


/**
 *  Partition a positive integer into its unique partitions.
 *
 *  @param n        the integer whose partitions are sought
 *  @param k        the partition size
 *
 *  @return a vector of the k-sized partitions
 * 
 *  @example Some examples of unique 3-way partitions are:
 *              {3, 0, 0} - {2, 1, 0} - {1, 1, 1}
 *              {4, 0, 0} - {3, 1, 0} - {2, 2, 0} - {2, 1, 1}
 */
std::vector<std::vector<size_t>> generateUniquePartitionsOf(const size_t n, const size_t k) {

    // The main algorithm starts from {n, 0, ..., 0} and moves a 1 from the right-most number (>1) to the left-most position that holds a value at least 2 smaller. If there are none such numbers left, the algorithm is finished.
    // Note that:
    //  - The largest number L in the partition will always be at the first position; the partition will always be sorted from largest to smallest value.
    //  - The right-most 1s can be ignored, since moving a 1 cannot create a new partition.
    //  - Moving a 1 to a position that holds a value L-1 or higher can be ignored, since it will always create a permutation of a previous partition.


    std::vector<size_t> partition(k, 0);
    partition.reserve(k);

    partition[0] = n;  // start with {n, 0, ..., 0}
    std::vector<std::vector<size_t>> unique_partitions {partition};
    while (true) {

        // Find the right-most position that holds a value larger than 1.
        auto subtraction_it = std::find_if(partition.rbegin(), partition.rend(), [](const size_t x) { return x > 1; });
        if (subtraction_it == partition.rend()) {
            break;
        }

        // Find the left-most position that is smaller than (value - 1).
        size_t value = *subtraction_it;
        auto addition_it = std::find_if(partition.begin(), partition.end(), [value](const size_t x) { return x < value - 1; });
        if (addition_it == partition.end()) {
            break;
        }

        // If there are such positions, proceed to move a 1.
        (*subtraction_it)--;
        (*addition_it)++;
        unique_partitions.push_back(partition);
    }

    return unique_partitions;
}


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
