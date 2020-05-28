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


#include "Basis/ScalarBasis/CartesianExponents.hpp"
#include "Mathematical/LinearCombination.hpp"
#include "Mathematical/ScalarFunction.hpp"


namespace GQCP {


/**
 *  A class representing a Cartesian Gaussian-type orbital (GTO), which is often referred to as a Cartesian 'primitive'
 *
 *  Mathematically speaking, a Cartesian GTO is a real-valued scalar function taking an Euclidean vector (3D-vector) as argument
 *
 *  Contracted GTOs can be expressed as linear combinations of GTOs: LinearCombination<CartesianGTO>
 */
class CartesianGTO: public ScalarFunction<double, double, 3> {
private:
    double gaussian_exponent;                // exponent of the exponential
    double N;                                // normalization factor
    CartesianExponents cartesian_exponents;  // exponents of (x-X), (y-Y), (z-Z)
    Vector<double, 3> m_center;              // center of the GTO (X, Y, Z)


public:
    // CONSTRUCTORS

    /**
     *  @param gaussian_exponent        the exponent of the exponential
     *  @param cartesian_exponents      the exponents of x, y and z
     *  @param center                   the center of the Cartesian GTO
     */
    CartesianGTO(const double gaussian_exponent, const CartesianExponents& cartesian_exponents, const Vector<double, 3>& center);

    /**
     *  Default constructor setting everything to zero
     */
    CartesianGTO();


    // OPERATORS

    /**
     *  @param r        the value at which the GTO should be evaluated
     *
     *  @return the value of the normalized GTO at the given position
     */
    double operator()(const Vector<double, 3>& r) const override;

    /**
     *  @param rhs      the right-hand side of the operator ==
     *
     *  @return whether two ScalarFunctions are equal
     */
    bool operator==(const CartesianGTO& rhs) const;


    // STATIC PUBLIC METHODS

    /**
     *  @param gaussian_exponent        the exponent of the GTO
     *  @param cartesian_exponent       the exponent of the Cartesian function x, y, z
     *
     *  @return one of the components of the total normalization factor
     */
    static double calculateNormalizationFactorComponent(const double gaussian_exponent, const size_t cartesian_exponent);

    /**
     *  @param gaussian_exponent        the exponent of the GTO
     *  @param cartesian_exponents      the exponents of the Cartesian functions x, y, z
     *
     *  @return the total normalization factor
     */
    static double calculateNormalizationFactor(const double gaussian_exponent, const CartesianExponents& cartesian_exponents);


    // PUBLIC METHODS

    /**
     *  @return the Cartesian exponents for this Cartesian GTO
     */
    const CartesianExponents& cartesianExponents() const { return this->cartesian_exponents; }

    /**
     *  @param direction        the Cartesian direction in which the derivative should be calculated
     *
     *  @return the derivative of this Cartesian GTO (with respect to the electronic coordinates) in the x-, y-, or z-direction
     */
    LinearCombination<double, CartesianGTO> calculateDerivative(const CartesianDirection direction) const;

    /**
     *  @return the gradient of this Cartesian GTO with respect to the electronic coordinates
     */
    Vector<LinearCombination<double, CartesianGTO>, 3> calculateGradient() const;

    /**
     *  @return the center of this Cartesian GTO
     */
    const Vector<double, 3>& center() const { return this->m_center; }

    /**
     *  @return the Gaussian exponent for this Cartesian GTO
     */
    double gaussianExponent() const { return this->gaussian_exponent; }

    /**
     *  @return the total normalization factor for this Cartesian GTO
     */
    double normalizationFactor() const { return this->N; }
};


}  // namespace GQCP
