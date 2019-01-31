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
#ifndef GQCP_NUMERICALDERIVATIVECALCULATOR_HPP
#define GQCP_NUMERICALDERIVATIVECALCULATOR_HPP


#include "BaseNumericalDerivativeCalculator.hpp"

#include "common.hpp"


namespace GQCP {

/**
 *  Class for calculating and storing derivatives and the associated function values for unary functions
 */
class NumericalDerivativeCalculator : public BaseNumericalDerivativeCalculator {
private:
    UnaryFunction uf;

public:
    // CONSTRUCTORS
    /**
     *  @param uf                    The Unary function to derive
     *  @param start                 starting parameter around which we derive
     *  @param step_size             step size for the numeric derivation approach
     */
    NumericalDerivativeCalculator (const UnaryFunction& uf, double start, double step_size);

    // PUBLIC METHODS
    /**
     *  Calculates all derivatives and function values up to the requested order
     *
     *  @param order                         highest requested order of derivative
     */
    void calculateDerivatives(size_t order);
};


}  // GQCP


#endif  // GQCP_NUMERICALDERIVATIVECALCULATOR_HPP
