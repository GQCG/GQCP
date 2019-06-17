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
#include "OrbitalOptimization/Localization/ERJacobiLocalizer.hpp"

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
void ERJacobiLocalizer::calculateJacobiCoefficients(const HamiltonianParameters<double>& ham_par, const size_t i, const size_t j) {

    const auto g = ham_par.get_g();  // two-electron integrals

    this->A = 0.25 * (2*g(i,i,j,j) + 4*g(i,j,i,j) - g(i,i,i,i) - g(j,j,j,j));
    this->B = -this->A;
    this->C = g(i,j,j,j) - g(i,i,i,j);
}


/**
 *  @param ham_par      the current Hamiltonian parameters
 *  @param i            the index of spatial orbital 1
 *  @param j            the index of spatial orbital 2
 *
 *  @return the angle for which the derivative of the scalar function after the Jacobi rotation is zero (and the second derivative is positive), using the current trigoniometric polynomial coefficients
 */
double ERJacobiLocalizer::calculateOptimalRotationAngle(const HamiltonianParameters<double>& ham_par, const size_t i, const size_t j) const {

    const double denominator = std::sqrt(std::pow(this->B, 2) + std::pow(this->C, 2));
    return 0.25 * std::atan2(this->C / denominator, this->B / denominator);  // atan(y/x) = std::atan2(y,x)
}


/**
 *  @param ham_par              the current Hamiltonian parameters
 *  @param jacobi_rot_par       the Jacobi rotation parameters
 * 
 *  @return the change in the value of the scalar function (i.e. the ER localization index) if the given Jacobi rotation parameters would be used to rotate the given Hamiltonian parameters
 */
double ERJacobiLocalizer::calculateScalarFunctionChange(const HamiltonianParameters<double>& ham_par, const JacobiRotationParameters& jacobi_rot_par) const {

    return - (this->A + std::sqrt(std::pow(this->B, 2) + std::pow(this->C, 2)));  // formulate as minimization problem
}



}  // namespace GQCP
