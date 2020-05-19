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


#include "Mathematical/Representation/Matrix.hpp"

#include <map>
#include <numeric>


namespace GQCP {


/**
 *  A slice of a matrix that only exists implicitly.
 * 
 *  If the full matrix is unnecessary to know, and only a certain slice of the matrix is of interest, this class implements operator() that can be used with the row and column indices of the full matrix.
 * 
 *  @tparam _Scalar             the scalar representation of one element of the encapsulating matrix
 */
template <typename _Scalar>
class ImplicitMatrixSlice {
public:
    using Scalar = _Scalar;


private:
    std::map<size_t, size_t> rows_implicit_to_dense;  // maps the row indices of the implicit matrix to the row indices of the dense representation of the slice
    std::map<size_t, size_t> cols_implicit_to_dense;  // maps the column indices of the implicit matrix to the column indices of the dense representation of the slice

    MatrixX<Scalar> M;  // the dense representation of the slice


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Initialize an ImplicitMatrixSlice's members.
     * 
     *  @param rows_implicit_to_dense           maps the row indices of the implicit matrix to the row indices of the dense representation of the slice
     *  @param cols_implicit_to_dense           maps the column indices of the implicit matrix to the column indices of the dense representation of the slice
     *  @param M                                the dense representation of the slice
     */
    ImplicitMatrixSlice(const std::map<size_t, size_t>& rows_implicit_to_dense, const std::map<size_t, size_t>& cols_implicit_to_dense, const MatrixX<Scalar>& M) :
        rows_implicit_to_dense {rows_implicit_to_dense},
        cols_implicit_to_dense {cols_implicit_to_dense},
        M {M} {

        // Check if the maps are consistent with the dense representation of the slice.
        if (this->rows_implicit_to_dense.size() != this->M.rows()) {
            throw std::invalid_argument("ImplicitMatrixSlice(const std::map<size_t, size_t>&, const std::map<size_t, size_t>&, const MatrixX<Scalar>&): The given dense representation of the slice does not have a compatible number of rows.");
        }

        if (this->cols_implicit_to_dense.size() != this->M.cols()) {
            throw std::invalid_argument("ImplicitMatrixSlice(const std::map<size_t, size_t>&, const std::map<size_t, size_t>&, const MatrixX<Scalar>&): The given dense representation of the slice does not have a compatible number of columns.");
        }
    }


    /**
     *  Initialize an ImplicitMatrixSlice's members, with a zero matrix for the dense representation of the slice.
     * 
     *  @param rows_implicit_to_dense           maps the row indices of the implicit matrix to the row indices of the dense representation of the slice
     *  @param cols_implicit_to_dense           maps the column indices of the implicit matrix to the column indices of the dense representation of the slice
     */
    ImplicitMatrixSlice(const std::map<size_t, size_t>& rows_implicit_to_dense, const std::map<size_t, size_t>& cols_implicit_to_dense) :
        ImplicitMatrixSlice(rows_implicit_to_dense, cols_implicit_to_dense,
                            MatrixX<Scalar>::Zero(rows_implicit_to_dense.size(), cols_implicit_to_dense.size())) {}


    /**
     *  A default constructor setting everything to zero.
     */
    ImplicitMatrixSlice() :
        // Use a named constructor for the default initialization.
        ImplicitMatrixSlice(ImplicitMatrixSlice<Scalar>::ZeroFromIndices({}, {})) {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Construct an ImplicitMatrixSlice from a dense matrix block and corresponding index ranges.
     * 
     *  @param row_start        the 0-based row index of the implicit matrix at which the block starts
     *  @param row_end          the 0-based row index of the implicit matrix at which the block ends (not included)
     *  @param col_start        the 0-based column index of the implicit matrix at which the block starts
     *  @param col_end          the 0-based column index of the implicit matrix at which the block ends (not included)
     *  @param M                the dense representation of the block
     * 
     *  @return an implicit matrix slice
     */
    static ImplicitMatrixSlice<Scalar> FromBlockRanges(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end, const MatrixX<Scalar>& M) {

        // Convert the ranges into allowed row and column indices of the implicit matrix, and then use another named constructor.
        std::vector<size_t> row_indices(row_end - row_start);
        std::iota(row_indices.begin(), row_indices.end(), row_start);

        std::vector<size_t> col_indices(col_end - col_start);
        std::iota(col_indices.begin(), col_indices.end(), col_start);

        return ImplicitMatrixSlice<Scalar>::FromIndices(row_indices, col_indices, M);
    }


    /**
     *  Create an implicit matrix slice through a dense representation of the slice.
     *
     *  @param row_indices              the row indices (in order) of the implicit matrix that the row indices of the dense representation of the slice correspond to
     *  @param col_indices              the column indices (in order) of the implicit matrix that the column indices of the dense representation of the slice correspond to
     *  @param M                        the dense representation of the slice
     * 
     *  @return an implicit matrix slice
     */
    static ImplicitMatrixSlice<Scalar> FromIndices(const std::vector<size_t>& row_indices, const std::vector<size_t>& col_indices, const MatrixX<Scalar>& M) {

        // Loop over all row and column indices, mapping them to the indices of the dense representation of the slice.
        std::map<size_t, size_t> rows_map;
        size_t dense_row_index = 0;  // row index in the dense representation of the slice
        for (const auto& implicit_row_index : row_indices) {

            rows_map[implicit_row_index] = dense_row_index;
            dense_row_index++;  // the dense indices are contiguous
        }

        std::map<size_t, size_t> cols_map;
        size_t dense_col_index = 0;  // column index in the dense representation of the slice
        for (const auto& implicit_col_index : col_indices) {

            cols_map[implicit_col_index] = dense_col_index;
            dense_col_index++;  // the dense indices are contiguous
        }

        return ImplicitMatrixSlice<Scalar>(rows_map, cols_map, M);
    }


    /**
     *  Construct a zero ImplicitMatrixSlice with given block ranges.
     * 
     *  @param row_start        the 0-based row index of the implicit matrix at which the block starts
     *  @param row_end          the 0-based row index of the implicit matrix at which the block ends (not included)
     *  @param col_start        the 0-based column index of the implicit matrix at which the block starts
     *  @param col_end          the 0-based column index of the implicit matrix at which the block ends (not included)
     *  @param M                the dense representation of the block
     * 
     *  @return a zero implicit matrix slice
     */
    static ImplicitMatrixSlice<Scalar> ZeroFromBlockRanges(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end) {

        // Create the zero dense representation of the slice and then use another named constructor.
        const auto rows = row_end - row_start;
        const auto cols = col_end - col_start;
        const MatrixX<Scalar> M = MatrixX<Scalar>::Zero(rows, cols);

        return ImplicitMatrixSlice<Scalar>::FromBlockRanges(row_start, row_end, col_start, col_end, M);
    }


    /**
     *  Create a zero-initialized implicit matrix slice from given allowed row and column indices.
     *
     *  @param row_indices              the row indices (in order) of the implicit matrix that the row indices of the dense representation of the slice correspond to
     *  @param col_indices              the column indices (in order) of the implicit matrix that the column indices of the dense representation of the slice correspond to
     * 
     *  @return an implicit matrix slice
     */
    static ImplicitMatrixSlice<Scalar> ZeroFromIndices(const std::vector<size_t>& row_indices, const std::vector<size_t>& col_indices) {

        // Zero-initialize a matrix with the number of rows and columns and use a different constructor.
        const MatrixX<Scalar> M = MatrixX<Scalar>::Zero(row_indices.size(), col_indices.size());

        return ImplicitMatrixSlice<Scalar>::FromIndices(row_indices, col_indices, M);
    }


    /*
     *  OPERATORS
     */

    /**
     *  Access an element of this implicit matrix slice.
     * 
     *  @param row                  the row number in the implicit encapsulating matrix
     *  @param col                  the column number in the implicit encapsulating matrix
     * 
     *  @return a read-only element of the implicit encapsulating matrix
     */
    Scalar operator()(const size_t row, const size_t col) const {

        // Map the implicit row and column indices to the row and column indices of the dense representation of this slice and access accordingly.
        const size_t i = this->denseIndexOfRow(row);
        const size_t j = this->denseIndexOfColumn(col);

        return this->M(i, j);
    }


    /**
     *  Access an element of this implicit matrix slice.
     * 
     *  @param row                  the row number in the implicit encapsulating matrix
     *  @param col                  the column number in the implicit encapsulating matrix
     * 
     *  @return a writable element of the impliti encapsulating matrix
     */
    Scalar& operator()(const size_t row, const size_t col) {

        // Map the implicit row and column indices to the row and column indices of the dense representation of this slice and access accordingly.
        const size_t i = this->denseIndexOfRow(row);
        const size_t j = this->denseIndexOfColumn(col);

        return this->M(i, j);
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return this as a (column-major) vector
     */
    VectorX<Scalar> asVector() const {
        return this->M.pairWiseReduce();
    }


    /**
     *  @return this as a matrix
     */
    const MatrixX<Scalar>& asMatrix() const {
        return this->M;
    }


    /**
     *  Convert an implicit column index to the column index in the dense representation of this slice.
     * 
     *  @param col                  the column number in the implicit encapsulating matrix
     * 
     *  @return the column index the dense representation of this slice.
     */
    size_t denseIndexOfColumn(const size_t col) const { return this->cols_implicit_to_dense.at(col); }


    /**
     *  Convert an implicit row index to the row index in the dense representation of this slice.
     * 
     *  @param row                  the row number in the implicit encapsulating matrix
     * 
     *  @return the row index the dense representation of this slice.
     */
    size_t denseIndexOfRow(const size_t row) const { return this->rows_implicit_to_dense.at(row); }
};


}  // namespace GQCP
