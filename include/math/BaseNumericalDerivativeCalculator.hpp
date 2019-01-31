// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef GQCP_BASENUMERICALDERIVATIVECALCULATOR_HPP
#define GQCP_BASENUMERICALDERIVATIVECALCULATOR_HPP


#include "common.hpp"


namespace GQCP {


/**
 *  Base class for calculating and storing derivatives and the associated function values for unary functions
 */
class BaseNumericalDerivativeCalculator {
protected:
    double start;
    double step_size;

    std::vector<double> function_values;
    std::vector<double> derivatives;

    // CONSTRUCTOR
    /**
     *  @param start                 starting parameter around which we derive
     *  @param step_size             step size for the numeric derivation approach
     */
    BaseNumericalDerivativeCalculator (double start, double step_size);

    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseNumericalDerivativeCalculator() = 0;

    // PROTECTED METHODS
    /**
     *  Numerically computes the derivative as : (-1)^(n+1) 1/(s^n) * sum^n_i=0 (-1)^(i+1) * (n)choose(i) * f(x + i*s)
     *
     *  @return the n-th order derivative
     */
    double calculateDerivative(size_t n) const ;

public:
    // GETTERS
    double get_derivative(size_t order) const { return this->derivatives[order];}
    double get_function_value(size_t order) const { return this->function_values[order];}

    // PUBLIC METHODS
    /**
     *  Calculates all derivatives and function values up to the requested order
     *
     *  @param order                         highest requested order of derivative
     */
    virtual void calculateDerivatives(size_t order) = 0;
};


}  // GQCP


#endif  // GQCP_BASENUMERICALDERIVATIVECALCULATOR_HPP
