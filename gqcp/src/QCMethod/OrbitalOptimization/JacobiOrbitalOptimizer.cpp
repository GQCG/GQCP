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

#include "QCMethod/OrbitalOptimization/JacobiOrbitalOptimizer.hpp"

#include <queue>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param dim                             the dimension of the orbital space that should be scanned. The valid orbital indices then are 0 ... dim (not included)
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
 */
JacobiOrbitalOptimizer::JacobiOrbitalOptimizer(const size_t dim, const double convergence_threshold, const size_t maximum_number_of_iterations) :
    dim {dim},
    BaseOrbitalOptimizer(convergence_threshold, maximum_number_of_iterations) {}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param sq_hamiltonian           The current Hamiltonian.
 * 
 *  @return The unitary transformation that will be used to rotate the current Hamiltonian into the next iteration.
 */
RTransformation<double> JacobiOrbitalOptimizer::calculateNewRotationMatrix(const RSQHamiltonian<double>& sq_hamiltonian) const {

    return RTransformation<double>::FromJacobi(this->optimal_jacobi_with_scalar.first, sq_hamiltonian.numberOfOrbitals());
}


/**
 *  @param sq_hamiltonian           the current Hamiltonian
 * 
 *  @return if the algorithm is considered to be converged
 */
bool JacobiOrbitalOptimizer::checkForConvergence(const RSQHamiltonian<double>& sq_hamiltonian) const {

    const double optimal_correction = optimal_jacobi_with_scalar.second;

    if (std::abs(optimal_correction) < this->convergence_threshold) {
        return true;
    } else {
        return false;
    }
}


/**
 *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
 */
void JacobiOrbitalOptimizer::prepareConvergenceChecking(const RSQHamiltonian<double>& sq_hamiltonian) {

    this->prepareJacobiSpecificConvergenceChecking(sq_hamiltonian);

    // Every Jacobi orbital optimizer should set a pair_type with the best Jacobi rotation.
    this->optimal_jacobi_with_scalar = this->calculateOptimalJacobiParameters(sq_hamiltonian);
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param sq_hamiltonian           the current Hamiltonian
 *  @param dim                      the dimension of the orbital space that should be scanned. The valid orbital indices then are 0 ... dim (not included)
 * 
 *  @return the optimal Jacobi rotation and the corresponding value for the scalar function that can be obtained when the Jacobi rotation would have taken place
 */
std::pair<JacobiRotation, double> JacobiOrbitalOptimizer::calculateOptimalJacobiParameters(const RSQHamiltonian<double>& sq_hamiltonian) {

    const auto& cmp = this->comparer();  // cmp: 'comparer'
    std::priority_queue<pair_type, std::vector<pair_type>, decltype(cmp)> queue {cmp};

    for (size_t q = 0; q < this->dim; q++) {
        for (size_t p = q + 1; p < this->dim; p++) {                  // loop over p>q
            this->calculateJacobiCoefficients(sq_hamiltonian, p, q);  // initialize the trigoniometric polynomial coefficients

            const double theta = this->calculateOptimalRotationAngle(sq_hamiltonian, p, q);
            const JacobiRotation jacobi_rotation {p, q, theta};

            const double E_change = this->calculateScalarFunctionChange(sq_hamiltonian, jacobi_rotation);

            queue.emplace(jacobi_rotation, E_change);  // construct a pair_type
        }
    }

    return queue.top();
}


/**
 *  @return the comparer functor that is used to compare two pair_types
 */
std::function<bool(const JacobiOrbitalOptimizer::pair_type&, const JacobiOrbitalOptimizer::pair_type&)> JacobiOrbitalOptimizer::comparer() const {

    return [this](const pair_type& lhs, const pair_type& rhs) {
        if (lhs.second < rhs.second) {
            return false;
        } else {
            return true;
        }
    };  // lambda function
}


}  // namespace GQCP
