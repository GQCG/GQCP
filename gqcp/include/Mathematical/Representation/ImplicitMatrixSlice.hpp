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

#include <unordered_map>


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
    size_t row_start;  // the 0-based row index of the full matrix at which the block starts, i.e. the start of the range of values that the first argument of operator() should accept
    size_t row_end;    // the 0-based row index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the first argument of operator() should accept

    size_t col_start;  // the 0-based column index of the full matrix at which the block starts, i.e. the start of the range of values that the second argument of operator() should accept
    size_t col_end;    // the 0-based column index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the second argument of operator() should accept

    std::unordered_map<size_t, size_t> rows_implicit_to_dense;  // maps the row indices of the implicit matrix to the row indices of the dense representation of the slice
    std::unordered_map<size_t, size_t> cols_implicit_to_dense;  // maps the column indices of the implicit matrix to the column indices of the dense representation of the slice


    std::vector<size_t> col_indices;  // the column indices (in order) of the implicit matrix that the column indices of the dense representation of the slice correspond to

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
    ImplicitMatrixSlice(const std::unordered_map<size_t, size_t>& rows_implicit_to_dense, const std::unordered_map<size_t, size_t>& cols_implicit_to_dense, const MatrixX<Scalar>& M) :
        rows_implicit_to_dense {rows_implicit_to_dense},
        cols_implicit_to_dense {cols_implicit_to_dense},
        M {M} {

        // Check if the maps are consistent with the dense representation of the slice.
        if (this->rows_implicit_to_dense.size() != this->M.rows()) {
            throw std::invalid_argument("ImplicitMatrixSlice(const std::unordered_map<size_t, size_t>&, const std::unordered_map<size_t, size_t>&, const MatrixX<Scalar>&): The given matrix does not have a compatible number of rows.");
        }

        if (this->cols_implicit_to_dense.size() != this->M.cols()) {
            throw std::invalid_argument("ImplicitMatrixSlice(const std::unordered_map<size_t, size_t>&, const std::unordered_map<size_t, size_t>&, const MatrixX<Scalar>&): The given matrix does not have a compatible number of columns.");
        }
    }


    /**
     *  Construct an ImplicitMatrixSlice with a given block matrix.
     * 
     *  @param row_start        the 0-based row index of the full matrix at which the block starts, i.e. the start of the range of values that the first argument of operator() should accept
     *  @param row_end          the 0-based row index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the first argument of operator() should accept
     *  @param col_start        the 0-based column index of the full matrix at which the block starts, i.e. the start of the range of values that the second argument of operator() should accept
     *  @param col_end          the 0-based column index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the second argument of operator() should accept
     *  @param M                the dense representation of the slice
     */
    ImplicitMatrixSlice(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end, const MatrixX<Scalar>& M) :
        row_start {row_start},
        row_end {row_end},
        col_start {col_start},
        col_end {col_end},
        M {M} {
        if (row_end - row_start != M.rows()) {
            throw std::invalid_argument("ImplicitMatrixSlice(const size_t, const size_t, const size_t, const size_t, const MatrixX<Scalar>&): The given matrix does not have a compatible number of rows.");
        }

        if (col_end - col_start != M.cols()) {
            throw std::invalid_argument("ImplicitMatrixSlice(const size_t, const size_t, const size_t, const size_t, const MatrixX<Scalar>&): The given matrix does not have a compatible number of columns.");
        }
    }

    /**
     *  Construct an ImplicitMatrixSlice with a zero block matrix.
     * 
     *  @param row_start        the 0-based row index of the full matrix at which the block starts, i.e. the start of the range of values that the first argument of operator() should accept
     *  @param row_end          the 0-based row index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the first argument of operator() should accept
     *  @param col_start        the 0-based column index of the full matrix at which the block starts, i.e. the start of the range of values that the second argument of operator() should accept
     *  @param col_end          the 0-based column index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the second argument of operator() should accept
     */
    ImplicitMatrixSlice(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end) :
        ImplicitMatrixSlice(row_start, row_end, col_start, col_end,
                            MatrixX<Scalar>::Zero(row_end - row_start, col_end - col_start)) {}


    /**
     *  A default constructor setting everything to zero.
     */
    ImplicitMatrixSlice() :
        ImplicitMatrixSlice(0, 0, 0, 0) {}


    // /**
    //  *  Create an implicit matrix slice through a dense representation of the slice.
    //  *
    //  *  @param row_indices              the row indices (in order) of the implicit matrix that the row indices of the dense representation of the slice correspond to
    //  *  @param col_indices              the column indices (in order) of the implicit matrix that the column indices of the dense representation of the slice correspond to
    //  *  @param M                        the dense representation of the slice
    //  */
    // ImplicitMatrixSlice(const std::vector<size_t>& row_indices, const std::vector<size_t>& col_indices, const MatrixX<Scalar>& M) :
    //     row_indices {row_indices},
    //     col_indices {col_indices},
    //     M {M} {}


    /*
     *  OPERATORS
     */

    /**
     *  Access an element of this implicit matrix slice.
     * 
     *  @param i                the row number in the implicit encapsulating matrix
     *  @param j                the column number in the implicit encapsulating matrix
     * 
     *  @return a read-only element of the implicit encapsulating matrix
     */
    Scalar operator()(const size_t i, const size_t j) const {

        // Map the implicit row and column indices to the row and column indices of the dense representation of this slice and access accordingly.
        const size_t row = this->denseIndexOfRow(i);
        const size_t col = this->denseIndexOfColumn(j);

        return this->M(row, col);
    }


    /**
     *  Access an element of this implicit matrix slice.
     * 
     *  @param i                the row number in the implicit encapsulating matrix
     *  @param j                the column number in the implicit encapsulating matrix
     * 
     *  @return a writable element of the impliti encapsulating matrix
     */
    Scalar& operator()(const size_t i, const size_t j) {

        // Map the implicit row and column indices to the row and column indices of the dense representation of this slice and access accordingly.
        const size_t row = this->denseIndexOfRow(i);
        const size_t col = this->denseIndexOfColumn(j);

        return this->M(row, col);
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
     *  @param j                the column number in the implicit encapsulating matrix
     * 
     *  @return the column index the dense representation of this slice.
     */
    size_t denseIndexOfColumn(const size_t j) const {

        return j - this->col_start;
    }


    /**
     *  Convert an implicit row index to the row index in the dense representation of this slice.
     * 
     *  @param i                the row number in the implicit encapsulating matrix
     * 
     *  @return the row index the dense representation of this slice.
     */
    size_t denseIndexOfRow(const size_t i) const {

        return i - this->row_start;
    }
};


}  // namespace GQCP
