// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include "Mathematical/Functions/CartesianDirection.hpp"
#include "Mathematical/Functions/CartesianExponents.hpp"
#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/**
 *  A class that can calculate overlap integrals over primitive Cartesian GTOs.
 */
class PrimitiveOverlapIntegralEngine {
public:
    // PUBLIC METHODS

    /**
     *  @param K                                the coordinate of the left primitive Cartesian GTO
     *  @param alpha                            the Gaussian exponent of the left primitive Cartesian GTO
     *  @param left_cartesian_exponents         the Cartesian exponents (x,y,z) of the left primitive Cartesian GTO
     *  @param L                                the coordinate of the right primitive Cartesian GTO
     *  @param beta                             the Gaussian exponent of the right primitive Cartesian GTO
     *  @param right_cartesian_exponents        the Cartesian exponent (x,y,z) of the right primitive Cartesian GTO
     */
    double calculate(const Vector<double, 3>& K, const double alpha, const CartesianExponents& left_cartesian_exponents, const Vector<double, 3>& L, const double beta, const CartesianExponents& right_cartesian_exponents);
};


}  // namespace GQCP
