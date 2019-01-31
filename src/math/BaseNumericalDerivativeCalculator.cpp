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
#include "math/BaseNumericalDerivativeCalculator.hpp"

#include "math.h"
#include "utilities/miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTOR
 */
BaseNumericalDerivativeCalculator::BaseNumericalDerivativeCalculator (double start, double step_size) : start(start), step_size(step_size) {}


/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
BaseNumericalDerivativeCalculator::~BaseNumericalDerivativeCalculator() {}


/*
 * PROTECTED METHODS
 */

/**
 *  Numerically computes the derivative as : (-1)^(n+1) 1/(s^n) * sum^n_i=0 (-1)^(i+1) * (n)choose(i) * f(x + i*s)
 *
 *  @return the n-th order derivative
 */
double BaseNumericalDerivativeCalculator::calculateDerivative(size_t n) const {
    double derivative = 0;
    double s = this->step_size;
    for (size_t i = 0; i < n+1; i++) {
        derivative += pow(-1,(i+1)) * binomialCoefficient(n,i) * this->get_function_value(i);
    }
    derivative /= pow(s, n);
    derivative *= pow(-1, n+1);
    return derivative;
}


}  // namespace GQCG








