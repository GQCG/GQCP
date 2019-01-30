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
#include "Localization/ERJacobiLocalizer.hpp"

#include <cmath>
#include <queue>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param N_P                              the number of electron pairs
 *  @param threshold                        the threshold for maximization on subsequent localization indices
 *  @param maximum_number_of_iterations     the maximum number of iterations for the localization algorithm
 */
ERJacobiLocalizer::ERJacobiLocalizer(size_t N_P, double threshold, size_t maximum_number_of_iterations) :
    BaseERLocalizer(N_P, threshold, maximum_number_of_iterations)
{}



/*
 *  PRIVATE METHODS
 */

/**
 *  Calculate the coefficients A, B, C for the Jacobi rotations
 *
 *  @param ham_par      the Hamiltonian parameters in an orthonormal basis
 *  @param i            the index of spatial orbital 1
 *  @param j            the index of spatial orbital 2
 */
void ERJacobiLocalizer::calculateJacobiCoefficients(const HamiltonianParameters& ham_par, size_t i, size_t j) {

    auto g = ham_par.get_g();  // two-electron integrals

    this->A = 0.25 * (2*g(i,i,j,j) + 4*g(i,j,i,j) - g(i,i,i,i) - g(j,j,j,j));
    this->B = -this->A;
    this->C = g(i,j,j,j) - g(i,i,i,j);
}


/**
 *  @param ham_par      the Hamiltonian parameters in an orthonormal basis
 *  @param i            the index of spatial orbital 1
 *  @param j            the index of spatial orbital 2
 *
 *  @return the angle which maximizes the Edmiston-Ruedenberg localization index for the orbitals i and j
 */
double ERJacobiLocalizer::calculateMaximizingRotationAngle(const HamiltonianParameters& ham_par, size_t i, size_t j) const {

    double denominator = std::sqrt(std::pow(this->B, 2) + std::pow(this->C, 2));
    return 0.25 * std::atan2(this->C / denominator, this->B / denominator);  // atan(y/x) = std::atan2(y,x)
}


/**
 *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) that contain the two-electron integrals upon which the Edmiston-Ruedenberg localization index is calculated
 *
 *  @return the maximal Edmiston-Ruedenberg for the current Jacobi coefficients A, B, C
 */
double ERJacobiLocalizer::calculateMaximalLocalizationIndex(const HamiltonianParameters& ham_par) const {

    double D = ham_par.calculateEdmistonRuedenbergLocalizationIndex(this->N_P);
    return D + this->A + std::sqrt(std::pow(this->B, 2) + std::pow(this->C, 2));
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Localize the Hamiltonian parameters by maximizing the Edmiston-Ruedenberg localization index, using the 'best' Jacobi rotation in every iteration step
 *
 *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) that should be localized
 */
void ERJacobiLocalizer::localize(HamiltonianParameters& ham_par) {

    while (!(this->is_converged)) {

        double D_old = ham_par.calculateEdmistonRuedenbergLocalizationIndex(this->N_P);

        // Find the Jacobi parameters (i,j,theta) that maximize the Edmiston-Ruedenberg localization index
        std::priority_queue<JacobiRotationLocalizationIndex> max_q;  // an ascending queue (on localization index) because we have implemented JacobiRotationLocalizationIndex::operator<

        for (size_t j = 0; j < this->N_P; j++) {
            for (size_t i = j+1; i < this->N_P; i++) {  // loop over i>j
                this->calculateJacobiCoefficients(ham_par, i, j);

                double theta = this->calculateMaximizingRotationAngle(ham_par, i, j);
                JacobiRotationParameters jacobi_rot_par {i, j, theta};
                double D_rotated = this->calculateMaximalLocalizationIndex(ham_par);

                max_q.emplace(JacobiRotationLocalizationIndex {jacobi_rot_par, D_rotated});
            }
        }

        // Rotate the Hamiltonian parameters using the maximizing Jacobi parameters
        auto maximizing_jacobi_parameters = max_q.top().jacobi_rotation_parameters;
        ham_par.rotate(maximizing_jacobi_parameters);


        // Check for convergence
        double D_new = max_q.top().index_after_rotation;  // the 'promised' localization index after rotation is equal to the one that would be calculated
        if (std::abs(D_new - D_old) < this->threshold) {
            this->is_converged = true;
        } else {
            this->iterations++;

            if (this->iterations == this->maximum_number_of_iterations) {
                throw std::runtime_error("The localization algorithm did not converge within the given maximum number of iterations.");
            }
        }

    }  // while not converged
}


}  // namespace GQCP
