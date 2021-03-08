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

#include "Basis/Integrals/McMurchieDavidsonCoefficient.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param K                One of the Cartesian components of the center of the left Cartesian GTO.
 *  @param a                The Gaussian exponent of the left Cartesian GTO
 *  @param L                One of the Cartesian components of the center of the right Cartesian GTO.
 *  @param b                The Gaussian exponent of the right Cartesian GTO.
 */
McMurchieDavidsonCoefficient::McMurchieDavidsonCoefficient(const double K, const double a, const double L, const double b) :
    K {K},
    L {L},
    a {a},
    b {b} {}


/*
 *  MARK: McMurchie-Davidson behavior
 */

/**
 *  @param i                The Cartesian exponent of the left Cartesian GTO.
 *  @param j                The Cartesian exponent of the right Cartesian GTO.
 *  @param t                The degree of the Hermite Gaussian.
 * 
 *  @return The value for the McMurchie-Davidson expansion coefficient E^{i,j}_t.
 */
double McMurchieDavidsonCoefficient::operator()(const int i, const int j, const int t) const {

    // Check if t is out of bounds: 0 <= t <= i+j.
    if ((t < 0) || (t > (i + j))) {
        return 0.0;
    }


    // Provide the base recurrence case.
    else if ((t == 0) && (i == 0) && (j == 0)) {
        return std::exp(-this->reducedExponent() * std::pow(this->distance(), 2));
    }

    // Do the recurrence for E^{i+1, j}_t.
    else if (j == 0) {
        return 1.0 / (2 * this->totalExponent()) * this->operator()(i - 1, j, t - 1) +
               (t + 1) * this->operator()(i - 1, j, t + 1) -
               this->b / this->totalExponent() * this->distance() * this->operator()(i - 1, j, t);
    }

    else {  // Do the recurrence for E^{i, j+1}_t.
        return 1.0 / (2 * this->totalExponent()) * this->operator()(i, j - 1, t - 1) +
               (t + 1) * this->operator()(i, j - 1, t + 1) +
               this->a / this->totalExponent() * this->distance() * this->operator()(i, j - 1, t);
    }
}


/*
 *  MARK: Gaussian overlap behavior
 */

/**
 *  @return The center of mass of the Gaussian overlap distribution.
 */
double McMurchieDavidsonCoefficient::centerOfMass() const {

    return (this->a * this->K + this->b * L) / this->totalExponent();
}


/**
 *  @return The reduced exponent of the Gaussian overlap distribution.
 */
double McMurchieDavidsonCoefficient::reducedExponent() const {

    return (this->a * this->b) / this->totalExponent();
}


}  // namespace GQCP
