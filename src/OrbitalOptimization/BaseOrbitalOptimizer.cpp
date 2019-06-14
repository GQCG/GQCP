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

/**
 *  @param oo_options               the options for orbital optimization
 */
BaseOrbitalOptimizer::BaseOrbitalOptimizer(std::shared_ptr<OrbitalOptimizationOptions> oo_options) :
    oo_options (std::move(oo_options))
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Optimize the Hamiltonian parameters by subsequently
 *      - checking for convergence (see checkForConvergence())
 *      - rotating the Hamiltonian parameters with a newly found rotation matrix (see calculateNewRotationMatrix())
 * 
 *  @param ham_par      the initial (guess for the) Hamiltonian parameters
 */
void BaseOrbitalOptimizer::optimize(HamiltonianParameters<double>& ham_par) {

    if (!ham_par.areOrbitalsOrthonormal()) {
        throw std::invalid_argument("BaseOrbitalOptimizer::optimize(HamiltonianParameters<double>&): The given Hamiltonian parameters do not belong to an orthonormal basis.");
    }

    size_t number_of_oo_iterations {0};
    while (this->prepareConvergenceChecking(ham_par), !this->checkForConvergence(ham_par)) {  // result of the comma operator is the second operand, if not converged
        const auto U = this->calculateNewRotationMatrix(ham_par);
        ham_par.rotate(U);

        number_of_oo_iterations++;
        if (number_of_oo_iterations > this->oo_options->maximumNumberOfIterations()) {
            throw std::runtime_error("BaseOrbitalOptimizer::optimize(HamiltonianParameters<double>&): The orbital optimization procedure did not converge in the given amount of iterations.");
        }
    }

    this->is_converged = true;
}


}  // namespace GQCP
