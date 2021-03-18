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

#include <cmath>


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

    return boost::math::hypergeometric_1F1(n + 0.5, n + 1.5, -x) / (2 * n + 1);
}


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
complex BoysFunction::operator()(const size_t n, const complex z) const {

    if ((std::abs(z) < 11.0) ||
        ((std::abs(z) < 20.0) && (z.real() >= 0.0))) {

        // Equation (25) in Shavitt 1963.
        complex series {};
        int i = 0;
        complex term {};
        do {
            term = std::pow(z, i) / std::tgamma(n + i + 1.5);
            series += term;
            i++;
        } while (std::abs(term) > 1.0e-15);

        return 0.5 * std::tgamma(n + 0.5) * std::exp(-z) * series;
    }

    else if ((std::abs(z) >= 20.0) && (z.real() >= 0.0)) {

        // Equation (32) in Shavitt 1963.
        complex series {};
        int i = 0;
        complex term {};
        do {
            term = 1.0 / (std::tgamma(n - i + 0.5) * std::pow(z, i));
            series += term;
            i++;
        } while (std::abs(term) > 1.0e-20);

        return std::tgamma(n + 0.5) / (2.0 * std::pow(z, n + 0.5)) -
               std::tgamma(n + 0.5) / (2.0 * z) * std::exp(-z) * series;
    }

    else {
        throw std::invalid_argument("BoysFunction::operator()(const size_t, const complex): The complex Boys function received an exceptional argument.");
    }
}


}  // namespace GQCP
