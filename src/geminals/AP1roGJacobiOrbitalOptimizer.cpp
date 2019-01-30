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
#include "geminals/AP1roGJacobiOrbitalOptimizer.hpp"

#include <cmath>
#include <queue>

#include "optimization/NewtonMinimizer.hpp"
#include "geminals/AP1roG.hpp"
#include "geminals/AP1roGPSESolver.hpp"

#include <boost/math/constants/constants.hpp>



namespace GQCP {


/*
 * CONSTRUCTORS
 */
/**
 *  @param N_P                                  the number of electron pairs
 *  @param ham_par                              Hamiltonian parameters in an orthonormal orbital basis
 *  @param oo_threshold                         the threshold on the convergence of the energy during the OO procedure
 *  @param maximum_number_of_oo_iterations      the maximum number of iterations during the OO procedure

 *  The initial guess for the geminal coefficients is zero
 */
AP1roGJacobiOrbitalOptimizer::AP1roGJacobiOrbitalOptimizer(size_t N_P, const HamiltonianParameters& ham_par, double oo_threshold, const size_t maximum_number_of_oo_iterations) :
    BaseAP1roGSolver(N_P, ham_par),
    oo_threshold (oo_threshold),
    maximum_number_of_oo_iterations (maximum_number_of_oo_iterations)
{}


/**
 *  @param molecule                             the molecule used for the AP1roG calculation
 *  @param ham_par                              Hamiltonian parameters in an orthonormal orbital basis
 *  @param oo_threshold                         the threshold on the convergence of the energy during the OO procedure
 *  @param maximum_number_of_oo_iterations      the maximum number of iterations during the OO procedure
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGJacobiOrbitalOptimizer::AP1roGJacobiOrbitalOptimizer(const Molecule& molecule, const HamiltonianParameters& ham_par, double oo_threshold, const size_t maximum_number_of_oo_iterations) :
    BaseAP1roGSolver(molecule, ham_par),
    oo_threshold (oo_threshold),
    maximum_number_of_oo_iterations (maximum_number_of_oo_iterations)
{}



/*
 *  PUBLIC METHODS
 */
/**
 *  Calculate the coefficients
 *      - A1, B1, C1            to be used in occupied-occupied rotations
 *      - A2, B2, C2, D2, E2    to be used in occupied-virtual rotations
 *      - A3, B3, C3            to be used in virtual-virtual rotations
 *
 *  @param p    the index of spatial orbital 1
 *  @param q    the index of spatial orbital 2
 *  @param G    the AP1roG geminal coefficients
 */
void AP1roGJacobiOrbitalOptimizer::calculateJacobiCoefficients(size_t p, size_t q, const AP1roGGeminalCoefficients& G) {

    OneElectronOperator h_SO = this->ham_par.get_h();
    TwoElectronOperator g_SO = this->ham_par.get_g();


    // Implementation of the Jacobi rotation coefficients with disjoint cases for p and q

    // Occupied-occupied rotations: if p <= N_P and q <= N_P for computers
    if ((p < this->N_P) && (q < this->N_P)) {

        this->A1 = 0.0;
        this->B1 = 0.0;
        this->C1 = 0.0;


        for (size_t b = this->N_P; b < this->K; b++) {
            this->A1 -= 0.5 * (g_SO(b,p,b,p) - g_SO(b,q,b,q)) * (G(p,b) - G(q,b));
            this->B1 += 0.5 * (g_SO(b,p,b,p) - g_SO(b,q,b,q)) * (G(p,b) - G(q,b));
            this->C1 += g_SO(b,p,b,q) * (G(q,b) - G(p,b));
        }
    }


    // Occupied-virtual rotations: if p > N_P and q <= N_P for computers
    else if ((p >= this->N_P) && (q < this->N_P)) {

        this->A2 = 0.0;
        this->B2 = 0.0;
        this->C2 = 0.0;
        this->D2 = 0.0;
        this->E2 = 0.0;


        this->A2 += h_SO(p,p) - h_SO(q,q) + 0.375 * (g_SO(p,p,p,p) + g_SO(q,q,q,q)) * (1 - G(q,p)) - 0.25 * g_SO(p,p,q,q) * (7 + G(q,p)) + 0.5 * g_SO(p,q,p,q) * (3 + G(q,p));
        this->B2 += h_SO(q,q) - h_SO(p,p) + 2 * g_SO(p,p,q,q) + 0.5 * (g_SO(p,p,p,p) + g_SO(q,q,q,q)) * (G(q,p) - 1) - g_SO(p,q,p,q) * (1 + G(q,p));
        this->C2 += 2 * h_SO(p,q) + (g_SO(p,p,p,q) - g_SO(p,q,q,q)) * (1 - G(q,p));
        this->D2 += 0.125 * (g_SO(p,p,p,p) + g_SO(q,q,q,q) - 2 * (g_SO(p,p,q,q) + 2 * g_SO(p,q,p,q))) * (1 - G(q,p));
        this->E2 += 0.5 * (g_SO(p,p,p,q) - g_SO(p,q,q,q)) * (G(q,p) - 1);

        for (size_t j = 0; j < this->N_P; j++) {
            this->A2 += 2 * (g_SO(j,j,p,p) - g_SO(j,j,q,q)) - 0.5 * (g_SO(j,p,j,p) - g_SO(j,q,j,q)) * (2 + G(j,p));
            this->B2 += 2 * (g_SO(j,j,q,q) - g_SO(j,j,p,p)) + 0.5 * (g_SO(j,p,j,p) - g_SO(j,q,j,q)) * (2 + G(j,p));
            this->C2 += 4 * g_SO(j,j,p,q) - g_SO(j,p,j,q) * (2 + G(j,p));
        }

        for (size_t b = this->N_P; b < this->K; b++) {
            this->A2 += 0.5 * (g_SO(b,p,b,p) - g_SO(b,q,b,q)) * G(q,b);
            this->B2 += 0.5 * (g_SO(b,q,b,q) - g_SO(b,p,b,p)) * G(q,b);
            this->C2 += g_SO(b,p,b,q) * G(q,b);
        }
    }


    // Virtual-virtual rotations: if p > N_P and q > N_P for computers
    else if ((p >= this->N_P) && (q >= this->N_P )) {

        this->A3 = 0.0;
        this->B3 = 0.0;
        this->C3 = 0.0;


        for (size_t j = 0; j < this->N_P; j++) {
            this->A3 -= 0.5 * (g_SO(j,p,j,p) - g_SO(j,q,j,q)) * (G(j,p) - G(j,q));
            this->B3 += 0.5 * (g_SO(j,p,j,p) - g_SO(j,q,j,q)) * (G(j,p) - G(j,q));
            this->C3 += g_SO(j,p,j,q) * (G(j,q) - G(j,p));
        }
    }


    else {  // this means that p <= q
        throw std::invalid_argument("The given p and q are invalid: p must be larger than q.");
    }

    this->are_calculated_jacobi_coefficients = true;
}


/**
 *  @param jacobi_rotation_parameters       the Jacobi parameters that specify a Jacobi rotation
 *  @param G                                the AP1roG geminal coefficients
 *
 *  @return the AP1roG energy after a Jacobi rotation using analytical formulas
 */
double AP1roGJacobiOrbitalOptimizer::calculateEnergyAfterJacobiRotation(const JacobiRotationParameters& jacobi_rotation_parameters, const AP1roGGeminalCoefficients& G) const {

    if (!this->are_calculated_jacobi_coefficients) {
        throw std::runtime_error("calculateEnergyAfterJacobiRotation: You haven't calculated the Jacobi coefficients yet. You should call AP1roG::calculateJacobiCoefficients before calling this function.");
    }
    // Note that this throw only warns the first time: it does not detect if there have been changes to the coefficients, for example after a recalculation


    size_t p = jacobi_rotation_parameters.get_p();
    size_t q = jacobi_rotation_parameters.get_q();
    double theta = jacobi_rotation_parameters.get_angle();


    // The formula I have derived is an energy CORRECTION due to the Jacobi rotation, so we initialize the rotated energy by the initial energy
    double E = calculateAP1roGEnergy(G, this->ham_par);

    // I've written everything in terms of cos(2 theta), sin(2 theta), cos(4 theta) and sin(4 theta)
    double c2 = std::cos(2 * theta);
    double s2 = std::sin(2 * theta);

    double c4 = std::cos(4 * theta);
    double s4 = std::sin(4 * theta);


    // Implementation of the Jacobi rotated energy with disjoint cases for p and q

    // Occupied-occupied rotations: if p <= N_P and q <= N_P for computers
    if ((p < this->N_P) && (q < this->N_P)) {
        return E + (this->A1 + this->B1 * c2 + this->C1 * s2);
    }

    // Occupied-virtual rotations: if p > N_P and q <= N_P for computers
    else if ((p >= this->N_P) && (q < this->N_P)) {
        return E + (this->A2 + this->B2 * c2 + this->C2 * s2 + this->D2 * c4 + this->E2 * s4);
    }

    // Virtual-virtual rotations: if p > N_P and q > N_P for computers
    else if ((p >= this->N_P) && (q >= this->N_P)) {

        return E + (this->A3 + this->B3 * c2 + this->C3 * s2);
    }

    else {  // this means that p <= q
        throw std::invalid_argument("The given p and q are invalid: p must be larger than q.");
    }
}


/**
 *  @param p    the index of spatial orbital 1
 *  @param q    the index of spatial orbital 2
 *  @param G    the AP1roG geminal coefficients
 *
 *  @return the angle for which the derivative of the energy after the Jacobi rotation is zero (and the second derivative is positive)
 */
double AP1roGJacobiOrbitalOptimizer::findOptimalRotationAngle(size_t p, size_t q, const AP1roGGeminalCoefficients& G) const {

    if (!this->are_calculated_jacobi_coefficients) {
        throw std::runtime_error("findOptimalRotationAngle: You haven't calculated the Jacobi coefficients yet. You should call AP1roG::calculateJacobiCoefficients before calling this function.");
    }
    // Note that this throw only warns the first time: it does not detect if there have been changes to the coefficients, for example after a recalculation


    // Implementation of the optimal Jacobi rotation angle with disjoint cases for p and q

    // Occupied-occupied rotations: if p <= N_P and q <= N_P for computers
    if ((p < this->N_P) && (q < this->N_P)) {
        double denominator = std::sqrt(std::pow(this->B1, 2) + std::pow(this->C1, 2));
        return 0.5 * std::atan2(-this->C1 / denominator, -this->B1 / denominator);  // std::atan2(y,x) = tan^-1(y/x)
    }


    // Occupied-virtual rotations: if p > N_P and q <= N_P for computers
    else if ((p >= this->N_P) && (q < this->N_P)) {

        std::priority_queue<JacobiRotationEnergy> min_q;  // an ascending queue (on energy) because we have implemented the 'reverse' JacobiParameters::operator<

        // Construct a lambda gradient function
        VectorFunction gradient_function = [this](const Eigen::VectorXd& x) {
            Eigen::VectorXd gradient_vec (1);
            gradient_vec << (-2*this->B2 * std::sin(2*x(0)) + 2*this->C2 * std::cos(2*x(0)) - 4*this->D2 * std::sin(4*x(0)) + 4*this->E2 * std::cos(4*x(0)));
            return gradient_vec;
        };

        // Construct a lambda Hessian function
        MatrixFunction hessian_function = [this](const Eigen::VectorXd& x) {
            Eigen::MatrixXd hessian_matrix (1, 1);
            hessian_matrix << (-4*this->B2 * std::cos(2*x(0)) - 2*this->C2 * std::sin(2*x(0)) - 16*this->D2 * std::cos(4*x(0)) - 16*this->E2 * std::sin(4*x(0)));
            return hessian_matrix;
        };


        // Use three initial guesses to get the minimum
        auto half_pi = boost::math::constants::half_pi<double>();
        double quarter_pi = std::pow(half_pi, 2);
        std::vector<double> theta_values {0.0, half_pi, quarter_pi};
        for (const auto& theta : theta_values) {
            Eigen::VectorXd theta_vec (1);  // we can't implicitly convert a float to an Eigen::VectorXd so we make it ourselves
            theta_vec << theta;

            NewtonMinimizer minimizer (theta_vec, gradient_function, hessian_function);
            minimizer.solve();
            double theta_min = minimizer.get_solution()(0);  // get inside the Eigen::VectorXd

            JacobiRotationParameters jacobi_rot_par {p, q, theta_min};

            double E_min = this->calculateEnergyAfterJacobiRotation(jacobi_rot_par, G);
            min_q.emplace(JacobiRotationEnergy {jacobi_rot_par, E_min});
        }  // for theta

        Eigen::VectorXd theta_min_vec (1);  // we can't implicitly convert a float to an Eigen::VectorXd so we make it ourselves
        theta_min_vec << min_q.top().jacobi_rotation_parameters.get_angle();

        assert(hessian_function(theta_min_vec)(0,0) > 0);  // the Hessian of the minimal value of the three must be positive, otherwise we're not in a minimum

        return min_q.top().jacobi_rotation_parameters.get_angle();
    }


    // Virtual-virtual rotations: if p > N_P and q > N_P for computers
    else if ((p >= this->N_P) && (q >= this->N_P )) {
        double denominator = std::sqrt(std::pow(this->B3, 2) + std::pow(this->C3, 2));
        return 0.5 * std::atan2(-this->C3 / denominator, -this->B3 / denominator);  // std::atan2(y,x) = tan^-1(y/x)
    }


    else {  // this means that p <= q
        throw std::invalid_argument("The given p and q are invalid: p must be larger than q.");
    }
}


/**
 *  Optimize the AP1roG energy by consequently
 *      - solving the AP1roG equations
 *      - finding the optimal Jacobi transformation (i.e. the one that yields the lowest energy)
 */
void AP1roGJacobiOrbitalOptimizer::solve() {

    // Solve the PSEs before starting
    AP1roGPSESolver initial_pse_solver (this->N_P, this->ham_par);
    initial_pse_solver.solve();
    auto G = initial_pse_solver.get_geminal_coefficients();
    double E_old = calculateAP1roGEnergy(G, this->ham_par);

    size_t iterations = 0;
    while (!(this->is_converged)) {

        // Find the Jacobi parameters (p,q,theta) that minimize the energy
        std::priority_queue<JacobiRotationEnergy> min_q;  // an ascending queue (on energy) because we have implemented the 'reverse' JacobiRotationEnergy::operator<

        for (size_t q = 0; q < this->K; q++) {
            for (size_t p = q+1; p < this->K; p++) {  // loop over p>q
                this->calculateJacobiCoefficients(p, q, G);

                double theta = this->findOptimalRotationAngle(p, q, G);
                JacobiRotationParameters jacobi_rot_par {p, q, theta};
                double E_rotated = this->calculateEnergyAfterJacobiRotation(jacobi_rot_par, G);

                min_q.emplace(JacobiRotationEnergy {jacobi_rot_par, E_rotated});
            }
        }  // loop over p>q

        // Check for if there were no Jacobi parameters that lower the energy
        double E_promised = min_q.top().energy_after_rotation;
        if (E_promised > E_old) {
            std::cerr << "Did not find a rotation that would lower the energy." << std::endl;
        }


        auto optimal_jacobi_parameters = min_q.top().jacobi_rotation_parameters;

        // Using the found Jacobi parameters, rotate the basis with the corresponding orthogonal Jacobi matrix
        this->ham_par.rotate(optimal_jacobi_parameters);


        // Solve the PSEs in the rotated spatial orbital basis
        AP1roGPSESolver pse_solver (this->N_P, this->ham_par, G);  // use the unrotated solution G as initial guess for the PSEs in the rotated basis
        pse_solver.solve();
        G = pse_solver.get_geminal_coefficients();
        double E = calculateAP1roGEnergy(G, this->ham_par);

        // Check for convergence
        if (std::abs(E - E_old) < this->oo_threshold) {
            this->is_converged = true;

            // Set the solution
            this->geminal_coefficients = G;
            this->electronic_energy = calculateAP1roGEnergy(this->geminal_coefficients, this->ham_par);
        } else {
            iterations++;
            E_old = E;  // copy the current energy to be able to check for energy convergence.

            if (iterations == this->maximum_number_of_oo_iterations) {
                throw std::runtime_error("The orbital optimization procedure did not converge.");
            }
        }
    }  // while not converged
}


}  // namespace GQCP
