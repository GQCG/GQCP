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


#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/**
 *  An implementation of the McMurchie-Davidson expansion coefficients through recurrence relations.
 */
class McMurchieDavidsonCoefficient {
private:
    double K;  // one of the components of the center of the left Cartesian GTO
    double L;  // one of the components of the center of the left Cartesian GTO

    double alpha;  // the Gaussian exponent of the left Cartesian GTO
    double beta;   // the Gaussian exponent of the right Cartesian GTO


public:
    // CONSTRUCTORS

    /**
     *  @param K                one of the components of the center of the left Cartesian GTO
     *  @param alpha            the Gaussian exponent of the left Cartesian GTO
     *  @param L                one of the components of the center of the left Cartesian GTO
     *  @param beta             the Gaussian exponent of the right Cartesian GTO
     */
    McMurchieDavidsonCoefficient(const double K, const double alpha, const double L, const double beta);


    // OPERATORS

    /**
     *  @param i            the Cartesian exponent of the left Cartesian GTO
     *  @param j            the Cartesian exponent of the right Cartesian GTO
     *  @param t            the degree of the Hermite Gaussian
     * 
     *  @return the value for the McMurchie-Davidson expansion coefficient E^{i,j}_t
     */
    double operator()(const int i, const int j, const int t) const;


    // PUBLIC METHODS

    /**
     *  @return the center of mass of the Gaussian overlap distribution
     */
    double centerOfMass() const;

    /**
     *  @return (one component of) the distance vector between the left and right Cartesian Gaussian (left - right)
     */
    double distance() const { return this->K - this->L; }

    /**
     *  @return the reduced exponent of the Gaussian overlap distribution
     */
    double reducedExponent() const;

    /**
     *  @return the total exponent of the Gaussian overlap distribution
     */
    double totalExponent() const { return alpha + beta; }
};


}  // namespace GQCP
