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
 *  @param dim                             the dimension of the orbital space that should be scanned. The valid orbital indices then are 0 ... dim (not included)
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
 */
JacobiOrbitalOptimizer::JacobiOrbitalOptimizer(const size_t dim, const double convergence_threshold, const size_t maximum_number_of_iterations) :
    dim (dim),
    BaseOrbitalOptimizer(convergence_threshold, maximum_number_of_iterations)
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
 */
void JacobiOrbitalOptimizer::prepareConvergenceChecking(const HamiltonianParameters<double>& ham_par) {

    this->prepareJacobiSpecificConvergenceChecking(ham_par);

    // Every Jacobi orbital optimizer should set a pair_type with the best Jacobi rotation parameters
    this->optimal_jacobi_with_scalar = this->calculateOptimalJacobiParameters(ham_par);
}


/**
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return if the algorithm is considered to be converged
 */
bool JacobiOrbitalOptimizer::checkForConvergence(const HamiltonianParameters<double>& ham_par) const {

    const double optimal_correction = optimal_jacobi_with_scalar.second;

    if (std::abs(optimal_correction) < this->convergence_threshold) {
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
SquareMatrix<double> JacobiOrbitalOptimizer::calculateNewRotationMatrix(const HamiltonianParameters<double>& ham_par) const {
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

    const auto& cmp = this->comparer();
    std::priority_queue<pair_type, std::vector<pair_type>, decltype(cmp)> queue (cmp);

    for (size_t q = 0; q < this->dim; q++) {
        for (size_t p = q+1; p < this->dim; p++) {  // loop over p>q
            this->calculateJacobiCoefficients(ham_par, p,q);  // initialize the trigoniometric polynomial coefficients

            const double theta = this->calculateOptimalRotationAngle(ham_par, p,q);
            const JacobiRotationParameters jacobi_rot_par (p, q, theta);

            const double E_change = this->calculateScalarFunctionChange(ham_par, jacobi_rot_par);

            queue.emplace(jacobi_rot_par, E_change);  // construct a pair_type
        }
    }

    return queue.top();
}


/**
 *  @return the comparer functor that is used to compare two pair_types
 */
std::function<bool (const JacobiOrbitalOptimizer::pair_type&, const JacobiOrbitalOptimizer::pair_type&)> JacobiOrbitalOptimizer::comparer() const {

    return [this] (const pair_type& lhs, const pair_type& rhs) {
        if (lhs.second < rhs.second) {
            return false;
        } else {
            return true;
        }
    };  // lambda function
}


}  // namespace GQCP
