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


#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  A class that represents a transformation matrix between two orbital bases.
 * 
 *  @tparam TransformationScalar            the scalar representation of one of the elements of the transformation matrix
 */
template <typename _TransformationScalar>
class TransformationMatrix : public SquareMatrix<_TransformationScalar> {
public:
    using TransformationScalar = _TransformationScalar;


public:
    /*
     *  CONSTRUCTORS
     */

    using SquareMatrix<TransformationScalar>::SquareMatrix;  // inherit base constructors


    /*
     *  PUBLIC METHODS
     */

    /**
     *  In-place 'transform' this transformation matrix such that the resulting transformation matrix describes this and the other transformation together
     */
    void transform(const TransformationMatrix<TransformationScalar>& T) {

        (*this) = (*this) * T;
    }
};


}  // namespace GQCP
