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
#pragma once


#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/**
 *  A class that represents a block of an encapsulating matrix. If the full matrix is unnecessary to know, and only a block of the matrix is of interest, this class implements operator() that can be used with the row and column indices of the full matrix. 
 * 
 *  @tparam _Scalar             the scalar representation of one element of the encapsulating matrix
 */
template <typename _Scalar>
class BlockMatrix {
public:
    using Scalar = _Scalar;


private:
    size_t row_start;  // the 0-based row index of the full matrix at which the block starts, i.e. the start of the range of values that the first argument of operator() should accept
    size_t row_end;  // the 0-based row index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the first argument of operator() should accept

    size_t col_start;  // the 0-based column index of the full matrix at which the block starts, i.e. the start of the range of values that the second argument of operator() should accept
    size_t col_end;  // the 0-based column index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the second argument of operator() should accept

    MatrixX<Scalar> M;  // the matrix representation of the block


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Constructor with a given block matrix
     * 
     *  @param row_start        the 0-based row index of the full matrix at which the block starts, i.e. the start of the range of values that the first argument of operator() should accept
     *  @param row_end          the 0-based row index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the first argument of operator() should accept
     *  @param col_start        the 0-based column index of the full matrix at which the block starts, i.e. the start of the range of values that the second argument of operator() should accept
     *  @param col_end          the 0-based column index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the second argument of operator() should accept
     *  @param M                the block matrix
     */
    BlockMatrix(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end, const MatrixX<Scalar>& M) :
        row_start (row_start),
        row_end (row_end),
        col_start (col_start),
        col_end (col_end),
        M (M)
    {
        if (row_end - row_start != M.rows()) {
            throw std::invalid_argument("BlockMatrix(const size_t, const size_t, const size_t, const size_t, const MatrixX<Scalar>&): The given matrix does not have a compatible number of rows.");
        }

        if (col_end - col_start != M.cols()) {
            throw std::invalid_argument("BlockMatrix(const size_t, const size_t, const size_t, const size_t, const MatrixX<Scalar>&): The given matrix does not have a compatible number of columns.");
        }
    }

    /**
     *  Constructor that initializes a zero block matrix
     * 
     *  @param row_start        the 0-based row index of the full matrix at which the block starts, i.e. the start of the range of values that the first argument of operator() should accept
     *  @param row_end          the 0-based row index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the first argument of operator() should accept
     *  @param col_start        the 0-based column index of the full matrix at which the block starts, i.e. the start of the range of values that the second argument of operator() should accept
     *  @param col_end          the 0-based column index of the full matrix at which the block ends (not included), i.e. the end (not included) of the range of values that the second argument of operator() should accept
     */
    BlockMatrix(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end) :
        BlockMatrix(row_start, row_end, col_start, col_end, MatrixX<Scalar>::Zero(row_end-row_start, col_end-col_start))
    {}



    /*
     *  OPERATORS
     */

    /**
     *  @param row      the row number in the encapsulating matrix
     *  @param col      the column number in the encapsulating matrix
     * 
     *  @return an element of the encapsulating matrix
     */
    Scalar operator()(const size_t row, const size_t col) const {

        const size_t row_block = row - this->row_start;  // the row index in the blocked matrix
        const size_t col_block = col - this->col_start;  // the columns index in the blocked matrix

        return this->M(row_block, col_block);
    }

    /**
     *  @param row      the row number in the encapsulating matrix
     *  @param col      the column number in the encapsulating matrix
     * 
     *  @return a writable element of the encapsulating matrix
     */
    Scalar& operator()(const size_t row, const size_t col) {

        const size_t row_block = row - this->row_start;  // the row index in the blocked matrix
        const size_t col_block = col - this->col_start;  // the columns index in the blocked matrix

        return this->M(row_block, col_block);
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
};


}  // namespace GQCP
