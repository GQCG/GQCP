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

#include "Mathematical/Functions/BoysFunction.hpp"

#include <boost/math/special_functions/hypergeometric_1F1.hpp>


namespace GQCP {


/*
 *  MARK: Function evaluation
 */

/**
 *  Calculate the value for the Boys function F_n(x).
 * 
 *  @param n        The degree of the Boys function.
 *  @param x        The argument for the Boys function.
 * 
 *  @return The value F_n(x).
 */
double BoysFunction::operator()(const size_t n, const double x) const {

    const auto n_ = static_cast<double>(n);
    return boost::math::hypergeometric_1F1(n_ + 0.5, n_ + 1.5, -x) / (2 * n_ + 1);
}


}  // namespace GQCP
