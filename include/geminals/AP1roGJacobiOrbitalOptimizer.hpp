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
#ifndef AP1roGJacobiOrbitalOptimizer_hpp
#define AP1roGJacobiOrbitalOptimizer_hpp


#include "geminals/BaseAP1roGSolver.hpp"


namespace GQCP {


/**
 *  A class that is used to find a minimum of the electronic energy for AP1roG wave functions with respect to rotations of the underlying spatial orbital basis.
 *
 *  By using analytical Jacobi rotations, and subsequently re-solving the AP1roG PSEs, a new orbital basis is found that (hopefully) results in a lower AP1roG energy.
 */
class AP1roGJacobiOrbitalOptimizer : public BaseAP1roGSolver {
private:
    // PRIVATE STRUCTS
    /**
     *  A struct that holds JacobiRotationParameters and an energy after rotation
     *
     *  Since operator< is implemented, "optimal parameters" can easily be found using a priority queue
     */
    struct JacobiRotationEnergy {

        JacobiRotationParameters jacobi_rotation_parameters;
        double energy_after_rotation;  // AP1roG energy after the rotation has taken place

        /**
         *  An operator< that can be used to achieve a minimum priority queue: the order of arguments is reversed
         *
         *  @param other    the other JacobiRotationEnergy parameters
         *
         *  @return if the energy of this is 'smaller' than other
         */
        bool operator< (const JacobiRotationEnergy& other) const {
            return this->energy_after_rotation > other.energy_after_rotation;
        }
    };


private:
    // PRIVATE PARAMETERS
    bool is_converged = false;
    double oo_threshold;  // the threshold used for OO: convergence is achieved when E_current - E_previous < oo_threshold
    size_t maximum_number_of_oo_iterations;

    bool are_calculated_jacobi_coefficients = false;
    double A1=0.0, B1=0.0, C1=0.0;  // Jacobi rotation coefficients for occupied-occupied rotations
    double A2=0.0, B2=0.0, C2=0.0, D2=0.0, E2=0.0;  // Jacobi rotation coefficients for occupied-virtual rotations
    double A3=0.0, B3=0.0, C3=0.0;  // Jacobi rotation coefficients for virtual-virtual rotations


public:
    // CONSTRUCTORS
    /**
     *  @param N_P                                  the number of electron pairs
     *  @param ham_par                              Hamiltonian parameters in an orthonormal orbital basis
     *  @param oo_threshold                         the threshold on the convergence of the energy during the OO procedure
     *  @param maximum_number_of_oo_iterations      the maximum number of iterations during the OO procedure

     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGJacobiOrbitalOptimizer(size_t N_P, const HamiltonianParameters& ham_par, double oo_threshold=1.0e-08, const size_t maximum_number_of_oo_iterations=128);

    /**
     *  @param molecule                             the molecule used for the OO-AP1roG calculation
     *  @param ham_par                              Hamiltonian parameters in an orthonormal orbital basis
     *  @param oo_threshold                         the threshold on the convergence of the energy during the OO procedure
     *  @param maximum_number_of_oo_iterations      the maximum number of iterations during the OO procedure
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGJacobiOrbitalOptimizer(const Molecule& molecule, const HamiltonianParameters& ham_par, double oo_threshold=1.0e-08, const size_t maximum_number_of_oo_iterations=128);


    // PUBLIC METHODS
    /**
     *  Calculate the coefficients
     *      - A1, B1, C1            to be used in occupied-occupied rotations
     *      - A2, B2, C2, D2, E2    to be used in occupied-virtual rotations
     *      - A3, B3, C3            to be used in virtual-virtual rotations
     *
     *  @param p    the index of spatial orbital 1
     *  @param q    the index of spatial orbital 2
     *  @param G    the AP1roG geminal coefficients
     */
    void calculateJacobiCoefficients(size_t p, size_t q, const AP1roGGeminalCoefficients& G);

    /**
     *  @param jacobi_rotation_parameters       the Jacobi parameters that specify a Jacobi rotation
     *  @param G                                the AP1roG geminal coefficients
     *
     *  @return the AP1roG energy after a Jacobi rotation using analytical formulas
     */
    double calculateEnergyAfterJacobiRotation(const JacobiRotationParameters& jacobi_rotation_parameters, const AP1roGGeminalCoefficients& G) const;

    /**
     *  @param p    the index of spatial orbital 1
     *  @param q    the index of spatial orbital 2
     *  @param G    the AP1roG geminal coefficients
     *
     *  @return the angle for which the derivative of the energy after the Jacobi rotation is zero (and the second derivative is positive)
     */
    double findOptimalRotationAngle(size_t p, size_t q, const AP1roGGeminalCoefficients& G) const;

    /**
     *  Optimize the AP1roG energy by consequently
     *      - solving the AP1roG equations
     *      - finding the optimal Jacobi transformation (i.e. the one that yields the lowest energy)
     */
    void solve() override;
};


}  // namespace GQCP


#endif /* AP1roGJacobiOrbitalOptimizer_hpp */
