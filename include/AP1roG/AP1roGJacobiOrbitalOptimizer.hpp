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
#ifndef AP1roGJacobiOrbitalOptimizer_hpp
#define AP1roGJacobiOrbitalOptimizer_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "AP1roG/AP1roG.hpp"
#include "JacobiRotationParameters.hpp"


namespace GQCP {


/**
 *  A class that is used to find a minimum of the electronic energy for AP1roG wave functions with respect to rotations of the underlying spatial orbital basis.
 *
 *  By using analytical Jacobi rotations, and subsequently re-solving the AP1roG PSEs, a new orbital basis is found that (hopefully) results in a lower AP1roG energy.
 */
class AP1roGJacobiOrbitalOptimizer {
private:
    // PRIVATE STRUCTS
    /**
     *  A struct that holds JacobiRotationParameters and a @param energy_after_rotation, operator< is implemented so "optimal parameters" can easily be
     *  found using a priority queue
     */
    struct JacobiRotationEnergy {

        GQCP::JacobiRotationParameters jacobi_rotation_parameters;
        double energy_after_rotation;  // AP1roG energy after the rotation has taken place

        /**
         *  An operator< that can be used to achieve a minimum priority queue: the order of arguments is reversed
         */
        bool operator< (const JacobiRotationEnergy& other) const {
            return this->energy_after_rotation > other.energy_after_rotation;
        }
    };


private:
    // PRIVATE PARAMETERS
    const size_t K;  // the number of special orbitals
    const size_t N_P;  // the number of electron pairs

    bool is_converged = false;
    const double oo_threshold;  // the threshold used for OO: convergence is achieved when E_current - E_previous < oo_threshold
    const size_t maximum_number_of_oo_iterations;


    GQCP::HamiltonianParameters ham_par;

    bool are_calculated_jacobi_coefficients = false;
    double A1=0.0, B1=0.0, C1=0.0;  // Jacobi rotation coefficients for occupied-occupied rotations
    double A2=0.0, B2=0.0, C2=0.0, D2=0.0, E2=0.0;  // Jacobi rotation coefficients for occupied-virtual rotations
    double A3=0.0, B3=0.0, C3=0.0;  // Jacobi rotation coefficients for virtual-virtual rotations

    GQCP::AP1roG solution;


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given number of electron pairs @param N_P, Hamiltonian parameters @param ham_par, a threshold for the orbital optimization @param oo_threshold and a @param maximum_number_of_oo_iterations
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGJacobiOrbitalOptimizer(size_t N_P, const GQCP::HamiltonianParameters& ham_par, double oo_threshold=1.0e-08, const size_t maximum_number_of_oo_iterations=128);

    /**
     *  Constructor based on a given @param molecule, Hamiltonian parameters @param ham_par, a threshold for the orbital optimization @param oo_threshold and a @param maximum_number_of_oo_iterations
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGJacobiOrbitalOptimizer(const GQCP::Molecule& molecule, const GQCP::HamiltonianParameters& ham_par, double oo_threshold=1.0e-08, const size_t maximum_number_of_oo_iterations=128);


    // GETTERS
    const GQCP::AP1roG& get_solution() const { return this->solution; }
    const GQCP::HamiltonianParameters& get_optimized_hamiltonian_parameters() const { return this->ham_par; }


    // PUBLIC METHODS
    /**
     *  Given the two indices of spatial orbitals @param p and @param q that will be Jacobi-rotated, and the geminal coefficients @param G,     calculate the coefficients (which are @members)
     *      - A1, B1, C1            to be used in occupied-occupied rotations
     *      - A2, B2, C2, D2, E2    to be used in occupied-virtual rotations
     *      - A3, B3, C3            to be used in virtual-virtual rotations
     */
    void calculateJacobiCoefficients(size_t p, size_t q, const GQCP::AP1roGGeminalCoefficients& G);

    /**
     *  Calculate the AP1roG energy given the geminal coefficients @param G after the application of a Jacobi rotation with the parameters @param jacobi_rotation_parameters
     */
    double calculateEnergyAfterJacobiRotation(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters, const GQCP::AP1roGGeminalCoefficients& G) const;

    /**
     *  Given a Jacobi pair @param p and @param q and the geminal coefficients @param G, @return the optimal rotation angle, i.e. the angle for which the derivative
     *  of the energy after the Jacobi rotation is zero (and the second derivative is positive).
     */
    double findOptimalRotationAngle(size_t p, size_t q, const GQCP::AP1roGGeminalCoefficients& G) const;

    /**
     *  Optimize the AP1roG energy by consequently
     *      - solving the AP1roG equations
     *      - finding the optimal Jacobi transformation (i.e. the one that yields the lowest energy)
     */
    void solve();
};


}  // namespace GQCP

#endif /* AP1roGJacobiOrbitalOptimizer_hpp */
