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


#include "Mathematical/Functions/LondonCartesianGTO.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param B            The homogeneous magnetic field appearing in the London modification.
 *  @param gto          The base Cartesian GTO.
 */
LondonCartesianGTO::LondonCartesianGTO(const HomogeneousMagneticField& B, const CartesianGTO& gto) :
    B {B},
    gto {gto} {}


/*
 *  MARK: Magnetic field
 */

/**
 *  @return The k-vector of the London plane wave, i.e. the value of the vector potential at the Gaussian center.
 */
Vector<double, 3> LondonCartesianGTO::kVector() const {

    const auto& K = this->gto.center();
    return this->magneticField().vectorPotentialAt(K);
}


/*
 *  MARK: Functional evaluation
 */

/**
 *  Evaluate the London prefactor at a given point in space.
 * 
 *  @param r            The point in space.
 * 
 *  @return The London plane wave phase factor at the given point.
 */
complex LondonCartesianGTO::phaseFactor(const Vector<double, 3>& r) const {

    using namespace GQCP::literals;
    return std::exp(-1.0_ii * this->kVector().dot(r));
}


/**
 *  Evaluate the London (Cartesian) GTO at a given point in space.
 * 
 *  @param r        The point in space.
 *
 *  @return The value of the London (Cartesian) GTO at the given point.
 */
complex LondonCartesianGTO::operator()(const Vector<double, 3>& r) const {

    return this->phaseFactor(r) * this->gto(r);
}


/**
 *  @param direction            the Cartesian direction in which the derivative should be calculated
 *
 *  @return the derivative of this London Cartesian GTO with respect to the position coordinate in the x-, y-, or z-direction
 */
 EvaluableLinearCombination<complex, LondonCartesianGTO> LondonCartesianGTO::calculatePositionDerivative(const CartesianDirection direction) const {

    using namespace GQCP::literals;

    // The formula consists of a part with the derivative of the plane wave, and a part with the derivative of the underlying GTO.
    // We start with the GTO derivative part, based on CartesianGTO::calculatePositionDerivative:

    // The formula is a sum of two parts: the derivative of the exponential and the derivative of the linear term (if applicable)

    // Derivative of the exponential (for eg component x): -2\alpha * x times original primitive 
    CartesianExponents exponential_derivative_exponents = this->gto.cartesianExponents();
    exponential_derivative_exponents.exponents[direction] += 1;
    CartesianGTO exponential_derivative_gto {this->gto.gaussianExponent(), exponential_derivative_exponents, this->gto.center()};
    // turn this into a london gto
    LondonCartesianGTO exponential_derivative_london_gto = {this->B, exponential_derivative_gto};
    // get coefficient for linear combination
    complex exponential_derivative_coefficient = -2 * this->gto.gaussianExponent();

    // add as first term to linear combination (one of three)
    EvaluableLinearCombination<complex, LondonCartesianGTO> lc {exponential_derivative_coefficient, exponential_derivative_london_gto};  // lc: linear combination


    // If the exponent in x, y or z is non-zero, there is an extra contribution of the linear term
    // i * 1/x with i the original exponent for x
    if (this->gto.cartesianExponents().value(direction) > 0) {

        CartesianExponents linear_derivative_exponents = this->gto.cartesianExponents();
        linear_derivative_exponents.exponents[direction] -= 1;

        CartesianGTO linear_derivative_gto(this->gto.gaussianExponent(), linear_derivative_exponents, this->gto.center());
        // again, turn into london gto
        LondonCartesianGTO linear_derivative_london_gto = {this->B, linear_derivative_gto};
        // get coefficient
        complex linear_derivative_coefficient = this->gto.cartesianExponents().value(direction);

        lc += EvaluableLinearCombination<complex, LondonCartesianGTO>(linear_derivative_coefficient, linear_derivative_london_gto);
    }

    // third term: simply original London GTO + i * k_x
    double k_component = this->kVector()[direction];
    complex plane_wave_derivative_coefficient = 1.0_ii * k_component;

    lc += EvaluableLinearCombination<complex, LondonCartesianGTO>(plane_wave_derivative_coefficient, *this);

    return lc;
}


/**
 *  @return the gradient of this London Cartesian GTO with respect to the position coordinate
 */
 Vector<EvaluableLinearCombination<complex, LondonCartesianGTO>, 3> LondonCartesianGTO::calculatePositionGradient() const {

    // Calculate the gradient for each of the Cartesian components.
    Vector<EvaluableLinearCombination<complex, LondonCartesianGTO>, 3> gradient;
    for (const auto& direction : {CartesianDirection::x, CartesianDirection::y, CartesianDirection::z}) {
        gradient(direction) = this->calculatePositionDerivative(direction);
    }

    return gradient;
}


}  // namespace GQCP
