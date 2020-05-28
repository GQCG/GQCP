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

#pragma once


#include "QCMethod/OrbitalOptimization/JacobiOrbitalOptimizer.hpp"
#include "QCModel/Geminals/AP1roGGeminalCoefficients.hpp"


namespace GQCP {


/**
 *  A class that is used to find a minimum of the electronic energy for AP1roG wave functions with respect to rotations of the underlying spatial orbital basis. By using analytical Jacobi rotations, and subsequently re-solving the AP1roG PSEs, a new orbital basis is found that results in a lower AP1roG energy
 * 
 *  Note that, since this algorithm does not use the PSE Lagrangian, a rotation can jump out of the AP1roG manifold, so this algorithm should be used with care
 */
class AP1roGJacobiOrbitalOptimizer: public JacobiOrbitalOptimizer {
private:
    size_t N_P;  // the number of electron pairs

    double A1 = 0.0, B1 = 0.0, C1 = 0.0;                      // coefficients for occupied-occupied rotations
    double A2 = 0.0, B2 = 0.0, C2 = 0.0, D2 = 0.0, E2 = 0.0;  // coefficients for occupied-virtual rotations
    double A3 = 0.0, B3 = 0.0, C3 = 0.0;                      // coefficients for virtual-virtual rotations

    double E;                     // the electronic energy
    AP1roGGeminalCoefficients G;  // the current geminal coefficients


public:
    // CONSTRUCTORS

    /**
     *  @param N_P                              the number of electron pairs
     *  @param K                                the number of spatial orbitals
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGJacobiOrbitalOptimizer(const size_t N_P, const size_t K, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128);

    /**
     *  @param G                                the initial geminal coefficients
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     */
    AP1roGJacobiOrbitalOptimizer(const AP1roGGeminalCoefficients& G, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  Calculate the trigoniometric polynomial coefficients for the given Jacobi rotation indices
     *
     *  @param p            the index of spatial orbital 1
     *  @param q            the index of spatial orbital 2
     */
    void calculateJacobiCoefficients(const SQHamiltonian<double>& sq_hamiltonian, const size_t p, const size_t q) override;

    /**
     *  @param sq_hamiltonian       the current Hamiltonian
     *  @param p                    the index of spatial orbital 1
     *  @param q                    the index of spatial orbital 2
     *
     *  @return the angle for which the derivative of the scalar function after the Jacobi rotation is zero (and the second derivative is positive), using the current trigoniometric polynomial coefficients
     */
    double calculateOptimalRotationAngle(const SQHamiltonian<double>& sq_hamiltonian, const size_t p, const size_t q) const override;

    /**
     *  @param sq_hamiltonian               the current Hamiltonian
     *  @param jacobi_rot_par               the Jacobi rotation parameters
     * 
     *  @return the change in the value of the scalar function (i.e. the AP1roG energy) if the given Jacobi rotation parameters would be used to rotate the given Hamiltonian
     */
    double calculateScalarFunctionChange(const SQHamiltonian<double>& sq_hamiltonian, const JacobiRotationParameters& jacobi_rot_par) const override;

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
     * 
     *  In the case of this uncoupled AP1roG Jacobi orbital optimizer, we should solve the AP1roG PSEs at the start at every iteration, using the current orbitals
     */
    void prepareJacobiSpecificConvergenceChecking(const SQHamiltonian<double>& sq_hamiltonian) override;


    // PUBLIC METHODS

    /**
     *  @return the electronic energy calculated by this orbital optimizer
     */
    double electronicEnergy() const { return this->E; }
};


}  // namespace GQCP
