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
#include "OrbitalOptimization/JacobiOrbitalOptimizer.hpp"

#include <queue>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param dim              the dimension of the orbital space that should be scanned. The valid orbital indices then are 0 ... dim (not included)
 *  @param oo_options       the options for orbital optimization
 */
JacobiOrbitalOptimizer::JacobiOrbitalOptimizer(const size_t dim, const OrbitalOptimizationOptions& oo_options) : 
    dim (dim),
    BaseOrbitalOptimizer(oo_options)
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return if the algorithm is considered to be converged
 */
bool JacobiOrbitalOptimizer::checkForConvergence(const HamiltonianParameters<double>& ham_par) {

    double old_value = this->calculateScalarFunction(ham_par);

    this->optimal_jacobi_with_scalar = this->calculateOptimalJacobiParameters(ham_par);
    double new_value = optimal_jacobi_with_scalar.second;

    if (std::abs(new_value - old_value) < this->oo_options.convergence_threshold) {
        return true;
    } else {
        return false;
    }
}


/**
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return a unitary matrix that will be used to rotate the current Hamiltonian parameters into the next iteration
 */
SquareMatrix<double> JacobiOrbitalOptimizer::calculateNewRotationMatrix(const HamiltonianParameters<double>& ham_par) {
    return SquareMatrix<double>::FromJacobi(this->optimal_jacobi_with_scalar.first, ham_par.get_K());
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param ham_par      the Hamiltonian parameters
 *  @param dim          the dimension of the orbital space that should be scanned. The valid orbital indices then are 0 ... dim (not included)
 * 
 *  @return the optimal Jacobi rotation parameters and the corresponding value for the scalar function that can be obtained when the Jacobi rotation would have taken place
 */
std::pair<JacobiRotationParameters, double> JacobiOrbitalOptimizer::calculateOptimalJacobiParameters(const HamiltonianParameters<double>& ham_par) {

    using pair_type = std::pair<JacobiRotationParameters, double>;

    auto cmp = [this] (const pair_type& lhs, const pair_type& rhs) {
        if (this->oo_options.should_minimize) {
            if (lhs.second < rhs.second) {
                return false;
            } else {
                return true;
            }
        }

        else {  // should maximize
            if (lhs.second < rhs.second) {
                return true;
            } else {
                return false;
            }
        }
    };

    std::priority_queue<pair_type, std::vector<pair_type>, decltype(cmp)> queue (cmp);


    for (size_t q = 0; q < this->dim; q++) {
        for (size_t p = q+1; p < this->dim; p++) {  // loop over p>q
            this->calculateJacobiCoefficients(ham_par, p,q);  // initialize the trigoniometric polynomial coefficients

            double theta = this->calculateOptimalRotationAngle(ham_par, p,q);
            const JacobiRotationParameters jacobi_rot_par (p, q, theta);

            double rotated_scalar_value = this->calculateScalarFunctionAfterJacobiRotation(ham_par, jacobi_rot_par);

            queue.emplace(jacobi_rot_par, rotated_scalar_value);  // constructs a pair_type
        }
    }

    return queue.top();
}




}  // namespace GQCP