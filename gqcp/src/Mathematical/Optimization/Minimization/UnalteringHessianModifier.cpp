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
#include "Mathematical/Optimization/Minimization/UnalteringHessianModifier.hpp"


namespace GQCP {


/**
 *  @param hessian      the current indefinite Hessian
 * 
 *  @return the given Hessian, i.e. do not alter the current hessian
 */
SquareMatrix<double> UnalteringHessianModifier::operator()(const SquareMatrix<double>& hessian) {
    return hessian;
}


}  // namespace GQCP
