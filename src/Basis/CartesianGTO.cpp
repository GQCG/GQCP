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
#include "Basis/CartesianGTO.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>



#include <iostream>

namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param alpha        the exponent of the exponential
 *  @param exponents    the exponents of x, y and z, in that order
 *  @param center       the center of the Cartesian GTO
 */
CartesianGTO::CartesianGTO(double alpha, const std::array<size_t, 3>& exponents, const Vector<double, 3>& center) :
    alpha (alpha),
    exponents (exponents),
    center (center)
{
    if (alpha < 0) {
        throw std::invalid_argument("CartesianGTO::CartesianGTO(double, std::array<size_t, 3>, Vector<double, 3>): the exponent must be larger than 0.");
    }


    this->N = this->calculateNormalizationFactor();
}


/**
 *  Default constructor setting everything to zero
 */
CartesianGTO::CartesianGTO() :
    CartesianGTO(0.0, std::array<size_t, 3> {0, 0, 0}, Vector<double, 3>::Zero())
{}



/*
 *  OPERATORS
 */

/**
 *  @param r        the value at which the GTO should be evaluated
 *
 *  @return the value of the GTO at the given position
 */
double CartesianGTO::operator()(const Vector<double, 3>& r) const {

    Vector<double, 3> delta_r = r - this->center;

    double value = this->N;
    value *= std::pow(delta_r.x(), this->exponents[0]);
    value *= std::pow(delta_r.y(), this->exponents[1]);
    value *= std::pow(delta_r.z(), this->exponents[2]);

    return value * std::exp(-this->alpha * delta_r.squaredNorm());
}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param alpha   the exponent of the GTO
 *  @param c       the power of the Cartesian function x, y, z
 *
 *  @return one of the components of the total normalization factor
 */
double CartesianGTO::calculateNormalizationFactorComponent(double alpha, size_t c) {

    double pi = boost::math::constants::pi<double>();

    double value = 1.0;
    value *= std::pow(2*alpha / pi, 0.25);

    if (c > 0) {  // if c==0, the following factor becomes 1
        auto df = boost::math::double_factorial<double>(2*static_cast<unsigned>(c) - 1);  // df: double factorial

        value *= std::pow(std::pow(4*alpha, c) / df, 0.5);
    }

    return value;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the normalization factor of the Cartesian GTO
 */
double CartesianGTO::calculateNormalizationFactor() const {
    double value = 1.0;

    // The formula is separable in its three components
    for (const auto& exponent : this->exponents) {
        value *= CartesianGTO::calculateNormalizationFactorComponent(this->alpha, exponent);
    }
    return value;
}


/**
 *  @param c        which component (x=0, y=1, z=2)
 *
 *  @return the derivative of this Cartesian GTO (with respect to the electronic coordinates) in the x-, y-, or z-direction
 */
LinearCombination<double, CartesianGTO> CartesianGTO::calculateDerivative(size_t c) const {

    // Derivative of the exponential
    std::array<size_t, 3> alpha_exponents = this->exponents;
    switch (c) {  // raise the exponents by one
        case 0:  // x-direction
            alpha_exponents[0] += 1;
            break;

        case 1:  // y-direction
            alpha_exponents[1] += 1;
            break;

        case 2:  // z-direction
            alpha_exponents[2] += 1;
            break;

        default:
            throw std::invalid_argument("CartesianGTO::calculateDerivative(size_t): the component can only be 0, 1, or 2");
            break;
    }

    CartesianGTO alpha_derivative (this->alpha, alpha_exponents, this->center);
    double alpha_coefficient = -2 * this->alpha;

    LinearCombination<double, CartesianGTO> lc (alpha_coefficient, alpha_derivative);  // lc: linear combination


    // If the exponent in x, y or z is non-zero, there is an extra contribution of the linear term
    if (this->exponents[c] > 0) {

        std::array<size_t, 3> linear_exponents = this->exponents;
        switch (c) {  // lower the exponents by one
            case 0:  // x-direction
                linear_exponents[0] -= 1;
                break;

            case 1:  // y-direction
                linear_exponents[1] -= 1;
                break;

            case 2:  // z-direction
                linear_exponents[2] -= 1;
                break;

            default:
                break;
        }

        CartesianGTO linear_derivative (this->alpha, linear_exponents, this->center);
        double linear_coefficient = this->exponents[c];

        lc += LinearCombination<double, CartesianGTO>(linear_coefficient, linear_derivative);
    }

    return lc;
}


/**
 *  @return the gradient of this Cartesian GTO with respect to the electronic coordinates
 */
Vector<LinearCombination<double, CartesianGTO>, 3> CartesianGTO::calculateGradient() const {

    Vector<LinearCombination<double, CartesianGTO>, 3> gradient;
    for (size_t c : {0, 1, 2}) {
        gradient(c) = this->calculateDerivative(c);
    }

    return gradient;
}


}  // namespace GQCP
