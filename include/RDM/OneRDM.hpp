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
#ifndef GQCP_ONERDM_HPP
#define GQCP_ONERDM_HPP


#include "math/SquareMatrix.hpp"


namespace GQCP {

/**
 *  A class that represents a 1-RDM
 *
 *  @tparam Scalar      the scalar type
 */
template <typename Scalar>
class OneRDM : public SquareMatrix<Scalar> {
public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  Default constructor
     */
    OneRDM() :
        SquareMatrix<Scalar>()
    {}


    /**
     *  @param matrix   the explicit matrix representation of the one-electron operator
     *
     *  Note that this should accept any Matrix<Scalar> (instead of SquareMatrix<Scalar>) because we want other Eigen return types to be accepted as well, like after a product of OneElectronOperators
     */
    explicit OneRDM(const Matrix<Scalar>& matrix) :
        SquareMatrix<Scalar>(matrix)
    {}
};


}  // namespace GQCP


#endif  // GQCP_ONERDM_HPP
