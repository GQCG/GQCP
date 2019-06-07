#include "OrbitalOptimization/NewtonOrbitalOptimizer.hpp"

#include "math/optimization/step.hpp"

#include <unsupported/Eigen/MatrixFunctions>


namespace GQCP {


/*
 *  PUBLIC OVERRIDDEN METHODS
 */ 

/**
 *  Determine if the algorithm has converged or not
 *  Specifically for the Newton-step based algorithms, this function
 *      - computes the gradient and checks its norm for convergence
 *      - if the gradient is zero, the Hessian is calculated and diagonalized and positive/negative definiteness is checked
 * 
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return if the algorithm is considered to be converged
 */
bool NewtonOrbitalOptimizer::checkForConvergence(const HamiltonianParameters<double>& ham_par) {

    std::cout << "Checking for convergence..." << std::endl;
    this->gradient = this->calculateGradientVector(ham_par);
    this->hessian = this->calculateHessianMatrix(ham_par);

    // Check for convergence on the norm
    std::cout << "norm of the gradient: " << this->gradient.norm() << std::endl;
    if (this->gradient.norm() < this->oo_options.convergence_threshold) {
        if (this->newtonStepIsWellDefined()) {
            return true;
        } else {
            return false;
        }
    }

    else {
        std::cout << "No convergence on the norm." << std::endl;
        return false;
    }
}


/**
 *  Produce a new rotation matrix by either
 *      - continuing in the direction of the largest (in absolute value) non-conforming eigenvalue (i.e. the smallest (negative) eigenvalue for minimization algorithms and the largest (positive) eigenvalue for maximization algorithms)
 *      - using the Newton step if it is well-defined
 * 
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return a unitary matrix that will be used to rotate the current Hamiltonian parameters into the next iteration
 */
SquareMatrix<double> NewtonOrbitalOptimizer::calculateNewRotationMatrix(const HamiltonianParameters<double>& ham_par) {

    // The general goal of this function is to:
    //      1) determine the free orbital rotation generators, using gradient and Hessian information
    //      2) construct the full orbital rotation generator matrix, by also including any redundant parameters
    //      3) calculate the unitary rotation matrix as the matrix exponential of the full orbital rotation generator matrix

    auto full_kappa_vector = this->calculateNewFullOrbitalGenerators(ham_par);  // should internally calculate the free orbital rotation generators

    auto full_kappa_matrix = GQCP::SquareMatrix<double>::FromStrictTriangle(full_kappa_vector);  // lower triangle only
    GQCP::SquareMatrix<double> full_kappa_matrix_transpose = full_kappa_matrix.transpose();
    full_kappa_matrix = full_kappa_matrix - full_kappa_matrix_transpose;  // add the antisymmetric component

    return (-full_kappa_matrix).exp();  // matrix exponential
}



/* 
 *  PUBLIC METHODS
 */ 

/**
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return the current orbital gradient as a vector
 */
VectorX<double> NewtonOrbitalOptimizer::calculateGradientVector(const HamiltonianParameters<double>& ham_par) const {
    return this->calculateGradientMatrix(ham_par).strictLowerTriangle();
}


/**
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return the current orbital Hessian as a matrix
 */
SquareMatrix<double> NewtonOrbitalOptimizer::calculateHessianMatrix(const HamiltonianParameters<double>& ham_par) const {
    return this->calculateHessianTensor(ham_par).pairWiseStrictReduce();
}


/**
 *  @return if a Newton step would be well-defined (i.e. the Hessian is positive definite for minimizations and negative definite for maximizations)
 */
bool NewtonOrbitalOptimizer::newtonStepIsWellDefined() const {

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hessian_diagonalizer (this->hessian);

    if (this->oo_options.should_minimize) {

        // Can only produce a well-defined descending Newton step if the Hessian is positive definite
        if (hessian_diagonalizer.eigenvalues()(0) < 0) {
            return false;
        } else {
            return true;
        }
    }
    
    else {  // should maximize

        // Can only produce a well-defined ascending Newton step if the Hessian is negative definite
        size_t dim = this->gradient.size();
        if (hessian_diagonalizer.eigenvalues()(dim - 1) > 0) {  // largest eigenvalue (they are sorted)
            return false;
        } else {
            return true;
        }
    }
}


/**
 *  If the Newton step is ill-defined, examine the Hessian and produce a new direction from it:
 *      - for minimization algorithms, this is the eigenvector that corresponds to the smallest (negative) eigenvalue of the Hessian
 *      - for maximization algorithms, this is the eigenvector that corresponds to the largest (positive) eigenvalue of the Hessian
 * 
 *  @return the new direction from the Hessian if the Newton step is ill-defined
 */
VectorX<double> NewtonOrbitalOptimizer::directionFromHessian() const {
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hessian_diagonalizer (this->hessian);

    if (this->oo_options.should_minimize) {
        if (hessian_diagonalizer.eigenvalues()(0) < 0) {
            return hessian_diagonalizer.eigenvectors().col(0);
        }
    }
    
    else {  // should maximize
        size_t dim = this->gradient.size();
        if (hessian_diagonalizer.eigenvalues()(dim - 1) > 0) {  // largest eigenvalue (they are sorted)
            return hessian_diagonalizer.eigenvectors().col(dim - 1);
        }
    }

    // We shouldn't ever reach these lines, unless this function was wrongly called
    throw std::invalid_argument("NewtonOrbitalOptimizer::directionFromHessian(): A direction from the Hessian was tried to be generated but the Newton step is well-defined.");
}


/**
 *  Use gradient and Hessian information to determine a new direction for the 'free' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
 * 
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return the new free orbital generators
 */
VectorX<double> NewtonOrbitalOptimizer::calculateNewFreeOrbitalGenerators(const HamiltonianParameters<double>& ham_par) const {

    // If the norm hasn't converged, continue in the Newton direction
    // Until we have completely read Nocedal & Wright, this is the way we're doing Newton optimization
    if (this->gradient.norm() > this->oo_options.convergence_threshold) {

        size_t dim = this->gradient.size();
        VectorFunction gradient_function = [this] (const VectorX<double>& x) { return this->gradient; };
        MatrixFunction hessian_function = [this] (const VectorX<double>& x) { return this->hessian; };

        return newtonStep(VectorX<double>::Zero(dim), gradient_function, hessian_function);  // with only the free parameters
    }

    else {  // the gradient has converged but the Hessian is indefinite
        // We're sure that if the program reaches this step, the Newton step is ill-defined
        std::cout << "Newton step is not well-defined, continuing in the Hessian direction..." << std::endl;
        return this->directionFromHessian();
    }
}


}  // namespace GQCP
