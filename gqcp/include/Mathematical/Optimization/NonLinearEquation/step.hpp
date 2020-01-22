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
#include "Utilities/typedefs.hpp"


namespace GQCP {


/**
 *  @param x        the current point
 *  @param f        a callable vector function
 *  @param J        the corresponding Jacobian function
 *
 *  @return the Newton step
 *      J(x) p = - f
 */
VectorX<double> newtonStep(const VectorX<double>& x, const VectorFunction<double>& f, const MatrixFunction<double>& J);


}  // namespace GQCP
