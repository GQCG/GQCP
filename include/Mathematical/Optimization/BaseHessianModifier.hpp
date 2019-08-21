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


namespace GQCP {


/**
 *  A base functor for Hessian modifiers
 */
class BaseHessianModifier {
public:
    // DESTRUCTOR
    virtual ~BaseHessianModifier() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  @param hessian      the current indefinite Hessian
     * 
     *  @return a modified Hessian that is made positive (for minimizers) or negative (for maximizers) definite
     */
    virtual SquareMatrix<double> operator()(const SquareMatrix<double>& hessian) = 0;
};


}  // namespace GQCP
