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
#include "OrbitalOptimization/BaseOrbitalOptimizer.hpp"


namespace GQCP {


/* 
 *  CONSTRUCTORS
 */

/*
*  @param convergence_threshold            the threshold used to check for convergence
*  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
*/
BaseOrbitalOptimizer::BaseOrbitalOptimizer(const double convergence_threshold, const size_t maximum_number_of_iterations) :
    convergence_threshold (convergence_threshold),
    maximum_number_of_iterations (maximum_number_of_iterations)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Optimize the Hamiltonian by subsequently
 *      - checking for convergence (see checkForConvergence())
 *      - rotating the Hamiltonian (and single-particle basis) with a newly found rotation matrix (see calculateNewRotationMatrix())
 * 
 *  @param sp_basis             the initial single-particle basis that contains the spinors to be optimized
 *  @param sq_hamiltonian       the initial (guess for the) Hamiltonian
 */
void BaseOrbitalOptimizer::optimize(SingleParticleBasis<double, GTOShell>& sp_basis, SQHamiltonian<double>& sq_hamiltonian) {

    if (!sp_basis.isOrthonormal()) {
        throw std::invalid_argument("BaseOrbitalOptimizer::optimize(SQHamiltonian<double>&): The given spinor basis is not orthonormal.");
    }

    while (this->prepareConvergenceChecking(sq_hamiltonian), !this->checkForConvergence(sq_hamiltonian)) {  // result of the comma operator is the second operand, so this expression effectively means "if not converged"
        const auto U = this->calculateNewRotationMatrix(sq_hamiltonian);
        sq_hamiltonian.rotate(U);

        this->number_of_iterations++;
        if (this->number_of_iterations > this->maximum_number_of_iterations) {
            throw std::runtime_error("BaseOrbitalOptimizer::optimize(SQHamiltonian<double>&): The orbital optimization procedure did not converge in the given amount of iterations.");
        }
    }

    this->is_converged = true;
}


}  // namespace GQCP
