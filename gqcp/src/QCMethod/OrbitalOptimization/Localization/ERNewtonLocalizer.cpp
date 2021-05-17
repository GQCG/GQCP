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

#include "QCMethod/OrbitalOptimization/Localization/ERNewtonLocalizer.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param orbital_space                    the occupied-virtual orbital space
 *  @param hessian_modifier                 the modifier functor that should be used when an indefinite Hessian is encountered
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
 */
ERNewtonLocalizer::ERNewtonLocalizer(const OrbitalSpace orbital_space, std::shared_ptr<BaseHessianModifier> hessian_modifier, const double convergence_threshold, const size_t maximum_number_of_iterations) :
    orbital_space {orbital_space},
    NewtonOrbitalOptimizer(hessian_modifier, convergence_threshold, maximum_number_of_iterations) {}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param sq_hamiltonian           the current Hamiltonian
 *
 *  @return the current orbital gradient of the Edmiston-Ruedenberg localization index as a matrix
 */
SquareMatrix<double> ERNewtonLocalizer::calculateGradientMatrix(const RSQHamiltonian<double>& sq_hamiltonian) const {

    const auto N_P = this->orbital_space.numberOfOrbitals(OccupationType::k_occupied);

    SquareMatrix<double> G = SquareMatrix<double>::Zero(N_P);

    for (const auto& i : this->orbital_space.indices(OccupationType::k_occupied)) {
        for (const auto& j : this->orbital_space.indices(OccupationType::k_occupied)) {
            G(i, j) = this->calculateGradientMatrixElement(sq_hamiltonian, i, j);
        }
    }

    return G;
}


/**
 *  @param sq_hamiltonian           the current Hamiltonian
 *
 *  @return the current orbital Hessian of the Edmiston-Ruedenberg localization index as a tensor
 */
SquareRankFourTensor<double> ERNewtonLocalizer::calculateHessianTensor(const RSQHamiltonian<double>& sq_hamiltonian) const {

    const auto N_P = this->orbital_space.numberOfOrbitals(OccupationType::k_occupied);
    SquareRankFourTensor<double> H {N_P};
    H.setZero();

    for (const auto& i : this->orbital_space.indices(OccupationType::k_occupied)) {
        for (const auto& j : this->orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& k : this->orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& l : this->orbital_space.indices(OccupationType::k_occupied)) {
                    H(i, j, k, l) = this->calculateHessianTensorElement(sq_hamiltonian, i, j, k, l);
                }
            }
        }
    }

    return H;
}


/**
 *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
 * 
 *  @param sq_hamiltonian           the current Hamiltonian
 * 
 *  @return the new full set orbital generators, including the redundant parameters
 */
ROrbitalRotationGenerators<double> ERNewtonLocalizer::calculateNewFullOrbitalGenerators(const RSQHamiltonian<double>& sq_hamiltonian) const {

    const auto kappa_free = this->calculateNewFreeOrbitalGenerators(sq_hamiltonian);  // only occupied-occupied
    const auto kappa_full = ROrbitalRotationGenerators<double>::FromOccOcc(kappa_free, sq_hamiltonian.numberOfOrbitals());

    return kappa_full;
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param sq_hamiltonian       the current Hamiltonian
 *  @param i                    the row of the gradient 'matrix'
 *  @param j                    the column of the gradient 'matrix'
 *
 *  @return the element (i,j) of the Edmiston-Ruedenberg localization index gradient
 */
double ERNewtonLocalizer::calculateGradientMatrixElement(const RSQHamiltonian<double>& sq_hamiltonian, const size_t i, const size_t j) const {

    const auto& g = sq_hamiltonian.twoElectron().parameters();

    return -4 * (g(j, i, i, i) - g(i, j, j, j));  // formulate as minimization problem
}


/**
 *  @param sq_hamiltonian       the current Hamiltonian
 *  @param i                    the first index of the Hessian 'tensor'
 *  @param j                    the second index of the Hessian 'tensor'
 *  @param k                    the third index of the Hessian 'tensor'
 *  @param l                    the fourth index of the Hessian 'tensor'
 *
 *  @return the element (i,j,k,l) of the Edmiston-Ruedenberg localization index Hessian
 */
double ERNewtonLocalizer::calculateHessianTensorElement(const RSQHamiltonian<double>& sq_hamiltonian, const size_t i, const size_t j, const size_t k, const size_t l) const {

    const auto& g = sq_hamiltonian.twoElectron().parameters();

    // KISS-implementation of the Hessian element for the Edmiston-Ruedenberg localization index
    double value = 0.0;
    if (i == k) {
        value += -2 * g(j, l, l, l) - 2 * g(l, j, j, j) + 8 * g(l, i, j, i) + 4 * g(l, j, i, i);
    }

    if (j == k) {
        value += 2 * g(i, l, l, l) + 2 * g(l, i, i, i) - 8 * g(l, j, i, j) - 4 * g(l, i, j, j);
    }

    if (i == l) {
        value += 2 * g(j, k, k, k) + 2 * g(k, j, j, j) - 8 * g(k, i, j, i) - 4 * g(k, j, i, i);
    }

    if (j == l) {
        value += -2 * g(i, k, k, k) - 2 * g(k, i, i, i) + 8 * g(k, j, i, j) + 4 * g(k, i, j, j);
    }

    return -value;  // formulate as minimization problem
}


}  // namespace GQCP
