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

#include "QCMethod/OrbitalOptimization/Localization/ERJacobiLocalizer.hpp"

#include <cmath>
#include <queue>


namespace GQCP {


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Calculate the trigoniometric polynomial coefficients for the given Jacobi rotation
 *
 *  @param i            the index of spatial orbital 1
 *  @param j            the index of spatial orbital 2
 */
void ERJacobiLocalizer::calculateJacobiCoefficients(const RSQHamiltonian<double>& sq_hamiltonian, const size_t i, const size_t j) {

    const auto& g = sq_hamiltonian.twoElectron().parameters();  // two-electron integrals

    this->A = 0.25 * (2 * g(i, i, j, j) + 4 * g(i, j, i, j) - g(i, i, i, i) - g(j, j, j, j));
    this->B = -this->A;
    this->C = g(i, j, j, j) - g(i, i, i, j);
}


/**
 *  @param sq_hamiltonian       the current Hamiltonian
 *  @param i                    the index of spatial orbital 1
 *  @param j                    the index of spatial orbital 2
 *
 *  @return the angle for which the derivative of the scalar function after the Jacobi rotation is zero (and the second derivative is positive), using the current trigoniometric polynomial coefficients
 */
double ERJacobiLocalizer::calculateOptimalRotationAngle(const RSQHamiltonian<double>& sq_hamiltonian, const size_t i, const size_t j) const {

    const double denominator = std::sqrt(std::pow(this->B, 2) + std::pow(this->C, 2));

    // If the denominator is almost zero, the Jacobi rotation is redundant: the corresponding angle of a 'non'-rotation is 0.0
    if (denominator < 1.0e-08) {
        return 0.0;
    }
    return 0.25 * std::atan2(this->C / denominator, this->B / denominator);  // atan(y/x) = std::atan2(y,x)
}


/**
 *  @param sq_hamiltonian               the current Hamiltonian
 *  @param jacobi_rot_par               the Jacobi rotation parameters
 * 
 *  @return the change in the value of the scalar function (i.e. the ER localization index) if the given Jacobi rotation parameters would be used to rotate the given Hamiltonian
 */
double ERJacobiLocalizer::calculateScalarFunctionChange(const RSQHamiltonian<double>& sq_hamiltonian, const JacobiRotationParameters& jacobi_rot_par) const {

    return -(this->A + std::sqrt(std::pow(this->B, 2) + std::pow(this->C, 2)));  // formulate as minimization problem
}


}  // namespace GQCP
