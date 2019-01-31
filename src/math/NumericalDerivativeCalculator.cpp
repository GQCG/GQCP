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
#include "math/NumericalDerivativeCalculator.hpp"


namespace GQCP {


/*
 *  CONSTRUCTOR
 */
NumericalDerivativeCalculator::NumericalDerivativeCalculator(const UnaryFunction& uf, double start, double step_size)
    : BaseNumericalDerivativeCalculator(start, step_size),
    uf(uf)
{

}


/*
 *  PUBLIC METHODS
 */

/**
 *  Calculates all derivatives and function values up to the requested order
 *
 *  @param order                         highest requested order of derivative
 */
void NumericalDerivativeCalculator::calculateDerivatives(size_t order) {
    for(size_t i = this->function_values.size(); i < order+1; i++) {
        this->function_values.push_back(uf(this->start + this->step_size*i));
    }

    for(size_t i = this->derivatives.size(); i < order+1; i++) {
        this->derivatives.push_back(this->calculateDerivative(i));
    }
}


}