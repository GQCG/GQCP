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
#ifndef ERJacobiLocalizer_hpp
#define ERJacobiLocalizer_hpp


#include "Localization/BaseERLocalizer.hpp"


namespace GQCP {


/**
 *  A class that localizes a set of orthonormal orbitals according to the maximization of the Edmiston-Ruedenberg localization index. A maximum is found using subsequent Jacobi rotations.
 */
class ERJacobiLocalizer : public BaseERLocalizer {
private:
    // PRIVATE MEMBERS
    bool are_calculated_jacobi_coefficients = false;
    double A=0.0, B=0.0, C=0.0;  // the Jacobi rotation coefficients


    // PRIVATE STRUCTS
    /**
     *  A struct that holds JacobiRotationParameters and a localization index
     *
     *  Since operator< is implemented, "optimal parameters" can easily be found using a priority queue
     */
    struct JacobiRotationLocalizationIndex {

        JacobiRotationParameters jacobi_rotation_parameters;
        double index_after_rotation;

        /**
         *  An operator< that can be used to achieve a minimum priority queue: the order of arguments is reversed
         *
         *  @param other    the other JacobiRotationLocalizationIndex parameters
         *
         *  @return if the localization index of this is smaller than other
         */
        bool operator< (const JacobiRotationLocalizationIndex& other) const {
            return this->index_after_rotation < other.index_after_rotation;
        }
    };


    // PRIVATE METHODS
    /**
     *  Calculate the coefficients A, B, C for the Jacobi rotations
     *
     *  @param ham_par      the Hamiltonian parameters in an orthonormal basis
     *  @param i            the index of spatial orbital 1
     *  @param j            the index of spatial orbital 2
     */
    void calculateJacobiCoefficients(const HamiltonianParameters& ham_par, size_t i, size_t j);

    /**
     *  @param ham_par      the Hamiltonian parameters in an orthonormal basis
     *  @param i            the index of spatial orbital 1
     *  @param j            the index of spatial orbital 2
     *
     *  @return the angle which maximizes the Edmiston-Ruedenberg localization index for the orbitals i and j
     */
    double calculateMaximizingRotationAngle(const HamiltonianParameters& ham_par, size_t i, size_t j) const;

    /**
     *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) that contain the two-electron integrals upon which the Edmiston-Ruedenberg localization index is calculated
     *
     *  @return the maximal Edmiston-Ruedenberg for the current Jacobi coefficients A, B, C
     */
    double calculateMaximalLocalizationIndex(const HamiltonianParameters& ham_par) const;


public:
    // CONSTRUCTORS
    /**
     *  @param N_P                              the number of electron pairs
     *  @param threshold                        the threshold for maximization on subsequent localization indices
     *  @param maximum_number_of_iterations     the maximum number of iterations for the localization algorithm
     */
    ERJacobiLocalizer(size_t N_P, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);


    // PUBLIC METHODS
    /**
     *  Localize the Hamiltonian parameters by maximizing the Edmiston-Ruedenberg localization index, using the 'best' Jacobi rotation in every iteration step
     *
     *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) that should be localized
     */
    void localize(HamiltonianParameters& ham_par) override;
};


}  // namespace GQCP


#endif /* ERJacobiLocalizer_hpp */
