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

#include "QCMethod/OrbitalOptimization/NewtonOrbitalOptimizer.hpp"

#include "Mathematical/Optimization/NonLinearEquation/step.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/*
 *  @param hessian_modifier                 the modifier functor that should be used when an indefinite Hessian is encountered
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
*/
NewtonOrbitalOptimizer::NewtonOrbitalOptimizer(std::shared_ptr<BaseHessianModifier> hessian_modifier, const double convergence_threshold, const size_t maximum_number_of_iterations) :
    hessian_modifier {hessian_modifier},
    BaseOrbitalOptimizer(convergence_threshold, maximum_number_of_iterations) {}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Produce a new rotation matrix by either
 *      - continuing in the direction of the largest (in absolute value) non-conforming eigenvalue (i.e. the smallest (negative) eigenvalue for minimization algorithms and the largest (positive) eigenvalue for maximization algorithms)
 *      - using the Newton step if it is well-defined
 * 
 *  @param sq_hamiltonian           the current Hamiltonian
 * 
 *  @return a unitary matrix that will be used to rotate the current Hamiltonian into the next iteration
 */
RTransformationMatrix<double> NewtonOrbitalOptimizer::calculateNewRotationMatrix(const SQHamiltonian<double>& sq_hamiltonian) const {

    // The general goal of this function is to:
    //      1) determine the free orbital rotation generators, using gradient and Hessian information
    //      2) determine the full orbital rotation generators, by also including any redundant parameters
    //      3) calculate the unitary rotation matrix from the full orbital rotation generators

    const auto full_kappa = this->calculateNewFullOrbitalGenerators(sq_hamiltonian);  // should internally calculate the free orbital rotation generators

    return full_kappa.calculateRotationMatrix();  // matrix exponential
}


/**
 *  Determine if the algorithm has converged or not
 *  Specifically for the Newton-step based algorithms, this function
 *      - computes the gradient and checks its norm for convergence
 *      - if the gradient is zero, the Hessian is calculated and positive definiteness is checked
 * 
 *  @param sq_hamiltonian           the current Hamiltonian
 * 
 *  @return if the algorithm is considered to be converged
 */
bool NewtonOrbitalOptimizer::checkForConvergence(const SQHamiltonian<double>& sq_hamiltonian) const {

    // Check for convergence on the norm
    if (this->gradient.norm() < this->convergence_threshold) {
        if (this->newtonStepIsWellDefined()) {  // needs this->hessian
            return true;
        } else {
            return false;
        }
    }

    else {
        return false;
    }
}


/**
 *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
 */
void NewtonOrbitalOptimizer::prepareConvergenceChecking(const SQHamiltonian<double>& sq_hamiltonian) {

    this->prepareOrbitalDerivativesCalculation(sq_hamiltonian);

    // All Newton-based orbital optimizers need to calculate a gradient and Hessian
    this->gradient = this->calculateGradientVector(sq_hamiltonian);
    this->hessian = this->calculateHessianMatrix(sq_hamiltonian);
}


/* 
 *  PUBLIC METHODS
 */

/**
 *  @param sq_hamiltonian           the current Hamiltonian
 * 
 *  @return the current orbital gradient as a vector
 */
VectorX<double> NewtonOrbitalOptimizer::calculateGradientVector(const SQHamiltonian<double>& sq_hamiltonian) const {
    return this->calculateGradientMatrix(sq_hamiltonian).pairWiseStrictReduced();
}


/**
 *  @param sq_hamiltonian           the current Hamiltonian
 * 
 *  @return the current orbital Hessian as a matrix
 */
SquareMatrix<double> NewtonOrbitalOptimizer::calculateHessianMatrix(const SQHamiltonian<double>& sq_hamiltonian) const {
    return this->calculateHessianTensor(sq_hamiltonian).pairWiseStrictReduced();
}


/**
 *  Use gradient and Hessian information to determine a new direction for the 'free' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
 * 
 *  @param sq_hamiltonian           the current Hamiltonian
 * 
 *  @return the new free orbital generators
 */
OrbitalRotationGenerators NewtonOrbitalOptimizer::calculateNewFreeOrbitalGenerators(const SQHamiltonian<double>& sq_hamiltonian) const {

    // If the norm hasn't converged, use the Newton step
    if (this->gradient.norm() > this->convergence_threshold) {

        const size_t dim = this->gradient.size();
        const VectorFunction<double> gradient_function = [this](const VectorX<double>& x) { return this->gradient; };

        auto modified_hessian = this->hessian;
        if (!this->newtonStepIsWellDefined()) {
            modified_hessian = this->hessian_modifier->operator()(this->hessian);
        }
        const MatrixFunction<double> hessian_function = [&modified_hessian](const VectorX<double>& x) { return modified_hessian; };

        return OrbitalRotationGenerators(newtonStep(VectorX<double>::Zero(dim), gradient_function, hessian_function));  // with only the free parameters
    }

    else {  // the gradient has converged but the Hessian is indefinite, so we have to 'push the algorithm over'
        // We're sure that if the program reaches this step, the Newton step is ill-defined
        return OrbitalRotationGenerators(this->directionFromIndefiniteHessian());
    }
}


/**
 *  If the Newton step is ill-defined, examine the Hessian and produce a new direction from it: the eigenvector that corresponds to the smallest (negative) eigenvalue of the Hessian
 * 
 *  @return the new direction from the Hessian if the Newton step is ill-defined
 */
VectorX<double> NewtonOrbitalOptimizer::directionFromIndefiniteHessian() const {

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hessian_diagonalizer {this->hessian};
    return hessian_diagonalizer.eigenvectors().col(0);
}


/**
 *  @return if a Newton step would be well-defined, i.e. the Hessian is positive definite
 */
bool NewtonOrbitalOptimizer::newtonStepIsWellDefined() const {

    // Can only produce a well-defined descending Newton step if the Hessian is positive definite
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hessian_diagonalizer {this->hessian};
    if (hessian_diagonalizer.eigenvalues()(0) < -1.0e-04) {  // not enough negative curvature to continue; can we change this to -this->convergence_threshold?
        return false;
    } else {
        return true;
    }
}


}  // namespace GQCP
