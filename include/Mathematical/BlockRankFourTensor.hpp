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
#ifndef GQCP_BLOCKRANKFOURTENSOR_HPP
#define GQCP_BLOCKRANKFOURTENSOR_HPP


#include "Mathematical/Matrix.hpp"
#include "Mathematical/Tensor.hpp"



namespace GQCP {


/**
 *  A class that represents a block of an encapsulating tensor
 * 
 *  @param _index1_start        the index of the first rank at which the block starts
 *  @param _index1_start        the index of the first rank at which the block ends (not included)
 * 
 *  @param _index2_start        the index of the second rank at which the block starts
 *  @param _index2_start        the index of the second rank at which the block starts (not included)
 * 
 *  @param _index3_start        the index of the third rank at which the block starts
 *  @param _index3_start        the index of the third rank at which the block starts (not included)
 * 
 *  @param _index4_start        the index of the fourth rank at which the block starts
 *  @param _index4_start        the index of the fourth rank at which the block starts (not included)
 */
template <typename _Scalar>
class BlockRankFourTensor {
public:
    using Scalar = _Scalar;


private:
    size_t index1_start;  // the index of the first rank at which the block starts
    size_t index1_end;  // the index of the first rank at which the block ends (not included)

    size_t index2_start;  // the index of the second rank at which the block starts
    size_t index2_end;  // the index of the second rank at which the block starts (not included)

    size_t index3_start;  // the index of the third rank at which the block starts
    size_t index3_end;  // the index of the third rank at which the block starts (not included)

    size_t index4_start;  // the index of the fourth rank at which the block starts
    size_t index4_end;  // the index of the fourth rank at which the block starts (not included)

    Tensor<Scalar, 4> T;  // the tensor representation of the block


public:
    /*
     *  CONSTRUCTOR
     */

    /**
     *  Constructor that initializes a zero block tensor
     * 
     *  @param index1_start         the index of the first rank at which the block starts
     *  @param index1_start         the index of the first rank at which the block ends (not included)
     * 
     *  @param index2_start         the index of the second rank at which the block starts
     *  @param index2_start         the index of the second rank at which the block starts (not included)
     * 
     *  @param index3_start         the index of the third rank at which the block starts
     *  @param index3_start         the index of the third rank at which the block starts (not included)
     * 
     *  @param index4_start         the index of the fourth rank at which the block starts
     *  @param index4_start         the index of the fourth rank at which the block starts (not included)
     */
    BlockRankFourTensor(const size_t index1_start, const size_t index1_end, const size_t index2_start, const size_t index2_end, const size_t index3_start, const size_t index3_end, const size_t index4_start, const size_t index4_end) :
        index1_start (index1_start),
        index1_end (index1_end),
        index2_start (index2_start),
        index2_end (index2_end),
        index3_start (index3_start),
        index3_end (index3_end),
        index4_start (index4_start),
        index4_end (index4_end),
        T (Tensor<Scalar, 4>(index1_end-index1_start, index2_end-index2_start, index3_end-index3_start, index4_end-index4_start))
    {
        T.setZero();
    }



    /*
     *  OPERATORS
     */

    /**
     *  @param index1           the first index of the encapsulating tensor
     *  @param index2           the second index of the encapsulating tensor
     *  @param index3           the third index of the encapsulating tensor
     *  @param index4           the fourth index of the encapsulating tensor
     * 
     *  @return an element of the encapsulating tensor
     */
    Scalar operator()(const size_t index1, const size_t index2, const size_t index3, const size_t index4) const {

        const size_t index1_block = index1 - this->index1_start;  // the first index in the blocked tensor
        const size_t index2_block = index2 - this->index2_start;  // the second index in the blocked tensor
        const size_t index3_block = index3 - this->index3_start;  // the third index in the blocked tensor
        const size_t index4_block = index4 - this->index4_start;  // the fourth index in the blocked tensor

        return this->T(index1_block,index2_block,index3_block,index4_block);
    }

    /**
     *  @param index1           the first index of the encapsulating tensor
     *  @param index2           the second index of the encapsulating tensor
     *  @param index3           the third index of the encapsulating tensor
     *  @param index4           the fourth index of the encapsulating tensor
     * 
     *  @return a writable element of the encapsulating tensor
     */
    Scalar& operator()(const size_t index1, const size_t index2, const size_t index3, const size_t index4) {

        const size_t index1_block = index1 - index1_start;  // the first index in the blocked tensor
        const size_t index2_block = index2 - index2_start;  // the second index in the blocked tensor
        const size_t index3_block = index3 - index3_start;  // the third index in the blocked tensor
        const size_t index4_block = index4 - index4_start;  // the fourth index in the blocked tensor

        return this->T(index1_block,index2_block,index3_block,index4_block);
    }



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return this as a (column-major) matrix
     */
    MatrixX<Scalar> asMatrix() const {
        return this->T.pairWiseReduce();
    }


    /**
     *  @return this as a tensor
     */
    const Tensor<double, 4>& asTensor() const {
        return this->T;
    }
};


}  // namespace GQCP



#endif  // GQCP_BLOCKRANKFOURTENSOR_HPP
