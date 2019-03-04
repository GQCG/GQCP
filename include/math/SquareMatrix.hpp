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
#ifndef SquareMatrix_hpp
#define SquareMatrix_hpp


#include "math/Matrix.hpp"


namespace GQCP {


/**
 *  A square extension of the Matrix class
 *
 *  @tparam Scalar      the scalar type
 */
template<typename Scalar>
class SquareMatrix : public MatrixX<Scalar> {
public:

    enum {
        Rows = MatrixX<Scalar>::RowsAtCompileTime,
        Cols = MatrixX<Scalar>::ColsAtCompileTime
    };


    /*
     *  CONSTRUCTORS
     */

    /**
     *  Default constructor
     */
    SquareMatrix() : MatrixX<Scalar>() {}


    /**
     *  A basic constructor that checks if the given matrix is square
     *
     *  @param matrix       the matrix that should be square
     */
    SquareMatrix(const MatrixX<Scalar>& matrix) :
        MatrixX<Scalar>(matrix)  // the compiler should call the move constructor here
    {
        // Check if the given matrix is square
        if (this->cols() != this->rows()) {
            throw std::invalid_argument("The given matrix is not square.");
        }
    }


    /*
     *  GETTERS
     */

    size_t get_dim() const {
        return this->cols();  // equals this->rows()
    }
};


}  // namespace GQCP


#endif /* SquareMatrix_hpp */
