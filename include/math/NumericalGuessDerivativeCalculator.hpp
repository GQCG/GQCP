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
#ifndef GQCP_NUMERICALGUESSDERIVATIVECALCULATOR_HPP
#define GQCP_NUMERICALGUESSDERIVATIVECALCULATOR_HPP


#include "BaseNumericalDerivativeCalculator.hpp"
#include "optimization/Eigenpair.hpp"

namespace GQCP {

using UnaryNumericEigenProblemFunction = std::function<Eigenpair (double x, const Eigen::VectorXd& guess)>;

/**
 *  Class for calculating and storing derivatives and the associated function values for unary functions
 */
class NumericalGuessDerivativeCalculator : public BaseNumericalDerivativeCalculator {
private:
    UnaryNumericEigenProblemFunction unepf;
    Eigen::VectorXd guess;
    std::vector<Eigen::VectorXd> eigenvectors;

public:
    // CONSTRUCTORS
    /**
     *  @param unepf                 The Unary function to derive
     *  @param start                 starting parameter around which we derive
     *  @param step_size             step size for the numeric derivation approach
     */
    NumericalGuessDerivativeCalculator (const UnaryNumericEigenProblemFunction& unepf, const Eigen::VectorXd& guess, double start, double step_size);


    // GETTERS
    const Eigen::VectorXd& get_eigenvector(size_t order) const { return this->eigenvectors[order]; }


    // PUBLIC METHODS
    /**
     *  Calculates all derivatives and function values up to the requested order
     *
     *  @param order                         highest requested order of derivative
     */
    void calculateDerivatives(size_t order);
};


}  // GQCP


#endif  // GQCP_NUMERICALGUESSDERIVATIVECALCULATOR_HPP
