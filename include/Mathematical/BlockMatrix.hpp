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
#ifndef GQCP_BLOCKMATRIX_HPP
#define GQCP_BLOCKMATRIX_HPP


#include "Mathematical/Matrix.hpp"


namespace GQCP {


/**
 *  A class that represents a block of an encapsulating matrix
 * 
 *  @tparam _Scalar             the scalar representation of one element of the encapsulating matrix
 */
template <typename _Scalar>
class BlockMatrix {
public:
    using Scalar = _Scalar;


private:
    size_t row_start;  // the 0-based row index at which the block starts
    size_t row_end;  // the 0-based row index at which the block ends (not included)

    size_t col_start;  // the 0-based column index at which the block starts
    size_t col_end;  // the 0-based column index at which the block ends (not included)

    MatrixX<Scalar> M;  // the matrix representation of the block


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Constructor that initializes a zero block matrix
     * 
     *  @param row_start        the 0-based row index at which the block starts
     *  @param row_end          the 0-based row index at which the block ends (not included)
     *  @param col_start        the 0-based column index at which the block starts
     *  @param col_end          the 0-based column index at which the block ends (not included)
     */
    BlockMatrix(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end) : 
        M (MatrixX<Scalar>::Zero(row_end-row_start, col_end-col_start))
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

        const size_t row_block = row - row_start;  // the row index in the blocked matrix
        const size_t col_block = col - col_start;  // the columns index in the blocked matrix

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



#endif  // GQCP_BLOCKMATRIX_HPP
