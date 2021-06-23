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


#include "Utilities/complex.hpp"

#include <cstddef>


namespace GQCP {


/**
 *  An implementation of the Boys function.
 */
class BoysFunction {
public:
    /*
     *  MARK: Function evaluation
     */

    /**
     *  Calculate the value for the real-valued Boys function F_n(x).
     * 
     *  @param n        The degree of the Boys function.
     *  @param x        The real-valued argument for the Boys function.
     * 
     *  @return The value F_n(x).
     */
    double operator()(const size_t n, const double x) const;

    /**
     *  Calculate the value for the complex-valued Boys function F_n(z).
     * 
     *  @param n        The degree of the Boys function.
     *  @param z        The real-valued argument for the Boys function.
     * 
     *  @return The value F_n(z).
     * 
     *  @note This method follows the implementation in Tellgren 2008, section II.C.
     */
    complex operator()(const size_t n, const complex z) const;
};


}  // namespace GQCP
