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
#include "Localization/ERNewtonLocalizer.hpp"

#include "utilities/linalg.hpp"
#include "math/optimization/step.hpp"


#include <unsupported/Eigen/MatrixFunctions>



namespace GQCP {


/*
 *  PRIVATE METHODS
 */
/**
 *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) containing the two-electron integrals
 *  @param i            the row of the gradient 'matrix'
 *  @param j            the column of the gradient 'matrix'
 *
 *  @return the element (i,j) of the Edmiston-Ruedenberg localization index gradient
 */
double ERNewtonLocalizer::calculateGradientElement(const HamiltonianParameters<double>& ham_par, size_t i, size_t j) const {

    auto g = ham_par.get_g();

    return 4 * (g(j,i,i,i) - g(i,j,j,j));
}


/**
 *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) containing the two-electron integrals
 *
 *  @return the gradient of the Edmiston-Ruedenberg localization index as a matrix
 */
Eigen::MatrixXd ERNewtonLocalizer::calculateGradient(const HamiltonianParameters<double>& ham_par) const {

    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(this->N_P, this->N_P);

    for (size_t i = 0; i < this->N_P; i++) {
        for (size_t j = 0; j < this->N_P; j++) {
            G(i,j) = this->calculateGradientElement(ham_par, i,j);
        }
    }

    return G;
}


/**
 *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) containing the two-electron integrals
 *  @param i            the first index of the Hessian 'tensor'
 *  @param j            the second index of the Hessian 'tensor'
 *  @param k            the third index of the Hessian 'tensor'
 *  @param l            the fourth index of the Hessian 'tensor'
 *
 *  @return the element (i,j,k,l) of the Edmiston-Ruedenberg localization index Hessian
 */
double ERNewtonLocalizer::calculateHessianElement(const HamiltonianParameters<double>& ham_par, size_t i, size_t j, size_t k, size_t l) const {

    auto g = ham_par.get_g();

    // KISS-implementation of the Hessian element for the Edmiston-Ruedenberg localization index
    double value = 0.0;
    if (i == k) {
        value += -2*g(j,l,l,l) - 2*g(l,j,j,j) + 8*g(l,i,j,i) + 4*g(l,j,i,i);
    }

    if (j == k) {
        value += 2*g(i,l,l,l) + 2*g(l,i,i,i) - 8*g(l,j,i,j) - 4*g(l,i,j,j);
    }

    if (i == l) {
        value += 2*g(j,k,k,k) + 2*g(k,j,j,j) - 8*g(k,i,j,i) - 4*g(k,j,i,i);
    }

    if (j == l) {
        value += -2*g(i,k,k,k) - 2*g(k,i,i,i) + 8*g(k,j,i,j) + 4*g(k,i,j,j);
    }

    return value;
}


/**
 *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) containing the two-electron integrals
 *
 *  @return the Hessian of the Edmiston-Ruedenberg localization index as a tensor
 */
Eigen::Tensor<double, 4> ERNewtonLocalizer::calculateHessian(const HamiltonianParameters<double>& ham_par) const {

    Eigen::Tensor<double, 4> H (this->N_P, this->N_P, this->N_P, this->N_P);
    H.setZero();

    for (size_t i = 0; i < this->N_P; i++) {
        for (size_t j = 0; j < this->N_P; j++) {
            for (size_t k = 0; k < this->N_P; k++) {
                for (size_t l = 0; l < this->N_P; l++) {
                    H(i,j,k,l) = this->calculateHessianElement(ham_par, i,j,k,l);
                }
            }
        }
    }

    return H;
}



/*
 *  CONSTRUCTORS
 */

/**
 *  @param N_P                              the number of electron pairs
 *  @param threshold                        the threshold for maximization on subsequent localization indices
 *  @param maximum_number_of_iterations     the maximum number of iterations for the localization algorithm
 */
ERNewtonLocalizer::ERNewtonLocalizer(size_t N_P, double threshold, size_t maximum_number_of_iterations) :
    BaseERLocalizer(N_P, threshold, maximum_number_of_iterations)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Localize the Hamiltonian parameters by maximizing the Edmiston-Ruedenberg localization index, using a Newton-based algorithm
 *
 *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) that should be localized
 */
void ERNewtonLocalizer::localize(HamiltonianParameters<double>& ham_par) {

    size_t dim = this->N_P * (this->N_P - 1) / 2 - 1;  // number of free occupied-occupied orbital rotation parameters

    while (!(this->is_converged)) {

        // Calculate the gradient and Hessian with only the free parameters, at kappa = 0
        Eigen::VectorXd gradient_vector = strictLowerTriangle(this->calculateGradient(ham_par));
        Eigen::MatrixXd hessian_matrix = strictLowerTriangle(this->calculateHessian(ham_par));

        // Perform a Newton-step to find orbital rotation parameters kappa
        VectorFunction gradient_function = [gradient_vector] (const Eigen::VectorXd& x) { return gradient_vector; };
        MatrixFunction hessian_function = [hessian_matrix] (const Eigen::VectorXd& x) { return hessian_matrix; };

        Eigen::VectorXd kappa_vector = newtonStep(Eigen::VectorXd::Zero(dim), gradient_function, hessian_function);  // with only the free parameters


        // If the calculated norm is zero, we have reached a critical point
        if (gradient_vector.norm() < this->threshold) {

            // If we have found a critical point, but we have a positive eigenvalue for the Hessian, continue in that direction
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hessian_solver (hessian_matrix);

            if (hessian_solver.eigenvalues()(dim - 1) > 0) {
                kappa_vector = hessian_solver.eigenvectors().col(dim - 1);
            }
            else {  // the Hessian is confirmed to be negative definite, so we have reached a maximum
                this->is_converged = true;
                break;
            }


        } else {
            this->iterations++;

            if (this->iterations == this->maximum_number_of_iterations) {
                throw std::runtime_error("The localization algorithm did not converge within the given maximum number of iterations.");
            }
        }


        // Change kappa from the occupied-occupied vector to the full matrix

        // The current kappa vector corresponds to the occupied-occupied rotation: the full kappa matrix contains o-o, o-v and v-v blocks
        Eigen::MatrixXd kappa_matrix_occupied = fillStrictLowerTriangle(kappa_vector);  // containing all parameters, so this is in anti-Hermitian (anti-symmetric) form
        Eigen::MatrixXd kappa_matrix_transpose_occupied = kappa_matrix_occupied.transpose();  // store the transpose in an auxiliary variable to avoid aliasing issues
        kappa_matrix_occupied -= kappa_matrix_transpose_occupied;  // fillStrictLowerTriangle only returns the lower triangle, so we must construct the anti-Hermitian (anti-symmetric) matrix

        Eigen::MatrixXd kappa_matrix = Eigen::MatrixXd::Zero(ham_par.get_K(), ham_par.get_K());
        kappa_matrix.topLeftCorner(this->N_P, this->N_P) = kappa_matrix_occupied;


        // Calculate the unitary rotation matrix that we can use to rotate the basis
        Eigen::MatrixXd U = (-kappa_matrix).exp();


        // Transform the integrals to the new orthonormal basis
        ham_par.Operator<double>::rotate(U);  // this checks if U is actually unitary
    }  // while not converged
}


}  // namespace GQCP
