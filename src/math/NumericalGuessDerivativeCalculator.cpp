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
#include "math/NumericalGuessDerivativeCalculator.hpp"


namespace GQCP {


/*
 *  CONSTRUCTOR
 */
NumericalGuessDerivativeCalculator::NumericalGuessDerivativeCalculator(const UnaryNumericEigenProblemFunction& unepf, const Eigen::VectorXd& guess, double start, double step_size)
    : BaseNumericalDerivativeCalculator(start, step_size),
      unepf(unepf),
      guess(guess)
{

}


/*
 *  PUBLIC METHODS
 */

/**
 *  Calculates all derivatives and function values and eigenvectors up to the requested order
 *
 *  @param order                         highest requested order of derivative
 */
void NumericalGuessDerivativeCalculator::calculateDerivatives(size_t order) {
    for(size_t i = function_values.size(); i < order+1; i++) {
        if (i == 0) {
            const Eigenpair& eigenpair = this->unepf(this->start, this->guess);
            this->function_values.push_back(eigenpair.get_eigenvalue());
            this->eigenvectors.push_back(eigenpair.get_eigenvector());
        } else {
            const Eigenpair& eigenpair = this->unepf(this->start + i*this->step_size, this->eigenvectors[i-1]);
            this->function_values.push_back(eigenpair.get_eigenvalue());
            this->eigenvectors.push_back(eigenpair.get_eigenvector());
        }
    }

    for(size_t i = derivatives.size(); i < order+1; i++) {
        this->derivatives.push_back(this->calculateDerivative(i));
    }
}


}