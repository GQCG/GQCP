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

#include "Mathematical/Functions/CartesianGTO.hpp"

#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <iostream>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Construct a Cartesian GTO from its members.
 * 
 *  @param gaussian_exponent        the exponent of the exponential
 *  @param cartesian_exponents      the exponents of x, y and z
 *  @param center                   the center of the Cartesian GTO
 */
CartesianGTO::CartesianGTO(const double gaussian_exponent, const CartesianExponents& cartesian_exponents, const Vector<double, 3>& center) :
    gaussian_exponent {gaussian_exponent},
    cartesian_exponents {cartesian_exponents},
    m_center {center} {

    if (gaussian_exponent < 0) {
        throw std::invalid_argument("CartesianGTO::CartesianGTO(const double, const CartesianExponents&, const Vector<double, 3>&): the Gaussian exponent must be larger than 0.");
    }
}


/**
 *  The default constructor setting everything to zero.
 */
CartesianGTO::CartesianGTO() :
    CartesianGTO(0.0, CartesianExponents(0, 0, 0), Vector<double, 3>::Zero()) {}


/*
 *  OPERATORS
 */

/**
 *  @param r        the value at which the GTO should be evaluated
 *
 *  @return the value of the unnormalized GTO at the given position
 */
double CartesianGTO::operator()(const Vector<double, 3>& r) const {

    Vector<double, 3> delta_r = r - this->m_center;

    double value = 1.0;
    value *= std::pow(delta_r.x(), this->cartesian_exponents.value(CartesianDirection::x));
    value *= std::pow(delta_r.y(), this->cartesian_exponents.value(CartesianDirection::y));
    value *= std::pow(delta_r.z(), this->cartesian_exponents.value(CartesianDirection::z));

    return value * std::exp(-this->gaussian_exponent * delta_r.squaredNorm());
}


/**
 *  @param other        the other CartesianGTO
 *
 *  @return if this CartesianGTO is equal to the other CartesianGTO
 */
bool CartesianGTO::operator==(const CartesianGTO& other) const {
    return (std::abs(this->gaussian_exponent - other.gaussian_exponent) < 1.0e-12) &&
           (this->cartesian_exponents == other.cartesian_exponents) &&
           (this->m_center.isApprox(other.m_center));
}


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param gaussian_exponent        the exponent of the GTO
 *  @param cartesian_exponent       the exponent of the Cartesian function (x, y, z)
 *
 *  @return one of the components (x, y, z) of the total normalization factor of a Cartesian GTO with the given exponents
 */
double CartesianGTO::calculateNormalizationFactorComponent(const double gaussian_exponent, const size_t cartesian_exponent) {

    double value = std::pow(2 * gaussian_exponent / boost::math::constants::pi<double>(), 0.25);

    // For a Cartesian exponent equal to zero, the factor becomes 1.
    if (cartesian_exponent > 0) {
        auto df = boost::math::double_factorial<double>(2 * static_cast<unsigned>(cartesian_exponent) - 1);  // df: double factorial
        value *= std::pow(std::pow(4 * gaussian_exponent, cartesian_exponent) / df, 0.5);
    }

    return value;
}


/**
 *  @param gaussian_exponent        the exponent of the GTO
 *  @param cartesian_exponents      the exponents of the Cartesian functions x, y, z
 *
 *  @return the total normalization factor of a Cartesian GTO with the given exponents
 */
double CartesianGTO::calculateNormalizationFactor(const double gaussian_exponent, const CartesianExponents& cartesian_exponents) {

    // The formula is separable in its three Cartesian components.
    double value = 1.0;
    for (const auto& exponent : cartesian_exponents.asArray()) {
        value *= CartesianGTO::calculateNormalizationFactorComponent(gaussian_exponent, exponent);
    }
    return value;
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param direction            the Cartesian direction in which the derivative should be calculated
 *
 *  @return the derivative of this Cartesian GTO with respect to the position coordinate in the x-, y-, or z-direction
 */
LinearCombination<double, CartesianGTO> CartesianGTO::calculatePositionDerivative(const CartesianDirection direction) const {

    // The formula is a sum of two parts: the derivative of the exponential and the derivative of the linear term (if applicable)

    // Derivative of the exponential
    CartesianExponents exponential_derivative_exponents = this->cartesian_exponents;
    exponential_derivative_exponents.exponents[direction] += 1;
    CartesianGTO exponential_derivative_gto {this->gaussian_exponent, exponential_derivative_exponents, this->m_center};
    double exponential_derivative_coefficient = -2 * this->gaussian_exponent;

    LinearCombination<double, CartesianGTO> lc {exponential_derivative_coefficient, exponential_derivative_gto};  // lc: linear combination


    // If the exponent in x, y or z is non-zero, there is an extra contribution of the linear term
    if (this->cartesian_exponents.value(direction) > 0) {

        CartesianExponents linear_derivative_exponents = this->cartesian_exponents;
        linear_derivative_exponents.exponents[direction] -= 1;

        CartesianGTO linear_derivative_gto(this->gaussian_exponent, linear_derivative_exponents, this->m_center);
        double linear_derivative_coefficient = this->cartesian_exponents.value(direction);

        lc += LinearCombination<double, CartesianGTO>(linear_derivative_coefficient, linear_derivative_gto);
    }

    return lc;
}


/**
 *  @return the gradient of this Cartesian GTO with respect to the position coordinate
 */
Vector<LinearCombination<double, CartesianGTO>, 3> CartesianGTO::calculatePositionGradient() const {

    // Calculate the gradient for each of the Cartesian components.
    Vector<LinearCombination<double, CartesianGTO>, 3> gradient;
    for (const auto& direction : {CartesianDirection::x, CartesianDirection::y, CartesianDirection::z}) {
        gradient(direction) = this->calculatePositionDerivative(direction);
    }

    return gradient;
}


/**
 *  @return a textual description of self
 */
std::string CartesianGTO::description() const {

    return (boost::format("CartesianGTO %s (%|.3f|)") % this->cartesian_exponents.description() % this->gaussian_exponent).str();
}


}  // namespace GQCP
