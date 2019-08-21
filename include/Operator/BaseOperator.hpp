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


#include "Mathematical/SquareMatrix.hpp"
#include "typedefs.hpp"


namespace GQCP {


/**
 *  An interface for second-quantized operators: they should implement the transformation formulas for their matrix representations in an orbital basis
 *
 *  CRTP is used for the static polymorphism, so the code will only compile if
 *      - DerivedOperator implements a suitable transform() method
 */
template <typename DerivedOperator>
class BaseOperator {
public:

    /**
     *  @return this as a DerivedOperator (done at compile time)
     */
    DerivedOperator& derived() { return static_cast<DerivedOperator&>(*this); }

    /**
     *  @return this as a const DerivedOperator (done at compile time)
     */
    const DerivedOperator& derived() const { return static_cast<DerivedOperator&>(*this); }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  In-place rotate the matrix representation of the operator
     *
     *  @param U     the unitary transformation (i.e. rotation) matrix, see basisTransform() for how the transformation matrix between the two bases should be represented
     */
    template<typename Scalar>
    void rotate(const SquareMatrix<Scalar>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("Operator::rotate(SquareMatrix<Scalar>): The given transformation matrix is not unitary.");
        }

        this->derived().basisTransform(U);
    }


    // Normally, I would have liked to put in the 'only-if-double' Jacobi rotation formula, but a function can't be and templated and virtual
};


}  // namespace GQCP
