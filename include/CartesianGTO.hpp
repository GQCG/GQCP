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
#ifndef CartesianGTO_hpp
#define CartesianGTO_hpp


#include "math/ScalarFunction.hpp"
#include "math/LinearCombination.hpp"


namespace GQCP {

/**
 *  A class representing a Cartesian Gaussian-type orbital (GTO), which is often referred to as a 'primitive'
 *
 *  Mathematically speaking, a Cartesian GTO is a real-valued scalar function taking an Euclidean vector (3D-vector) as argument
 *
 *  Contracted GTOs can be expressed as linear combinations of GTOs: LinearCombination<CartesianGTO>
 */
class CartesianGTO : public ScalarFunction<double, double, 3> {
public:
    double alpha;  // exponent of the exponential
    double N;  // normalization factor
    std::array<size_t, 3> exponents;  // exponents of (x-X), (y-Y), (z-Z)
    Vector<double, 3> center;  // center of the GTO (X, Y, Z)


public:
    // CONSTRUCTORS
    /**
     *  @param alpha        the exponent of the exponential
     *  @param exponents    the exponents of x, y and z, in that order
     *  @param center       the center of the Cartesian GTO
     */
    CartesianGTO(double alpha, const std::array<size_t, 3>& exponents, const Vector<double, 3>& center);

    /**
     *  Default constructor setting everything to zero
     */
    CartesianGTO();


    // GETTERS
    double get_exponent() const { return this->alpha; }
    const std::array<size_t, 3>& get_exponents() const { return this->exponents; }
    const Vector<double, 3>& get_center() const { return this->center; }


    // OPERATORS
    /**
     *  @param r        the value at which the GTO should be evaluated
     *
     *  @return the value of the GTO at the given position
     */
    double operator()(const Vector<double, 3>& r) const override;


    // STATIC PUBLIC METHODS
    /**
     *  @param alpha   the exponent of the GTO
     *  @param c       the power of the Cartesian function x, y, z
     *
     *  @return one of the components of the total normalization factor
     */
    static double calculateNormalizationFactorComponent(double alpha, size_t c);


    // PUBLIC METHODS
    /**
     *  @return the total normalization factor of the Cartesian GTO
     */
    double calculateNormalizationFactor() const;

    /**
     *  @param c        which component (x=0, y=1, z=2)
     *
     *  @return the derivative of this Cartesian GTO (with respect to the electronic coordinates) in the x-, y-, or z-direction
     */
    LinearCombination<double, CartesianGTO> calculateDerivative(size_t c) const;

    /**
     *  @return the gradient of this Cartesian GTO with respect to the electronic coordinates
     */
    Vector<LinearCombination<double, CartesianGTO>, 3> calculateGradient() const;
};


}  // namespace GQCP


#endif  /* CartesianGTO_hpp */
