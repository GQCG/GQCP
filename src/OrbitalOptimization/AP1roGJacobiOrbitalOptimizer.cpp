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
#include "OrbitalOptimization/AP1roGJacobiOrbitalOptimizer.hpp"

#include "Geminals/AP1roG.hpp"
#include "Geminals/AP1roGPSESolver.hpp"
#include "Mathematical/Optimization/NewtonMinimizer.hpp"

#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <queue>


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  @param N_P                              the number of electron pairs
 *  @param K                                the number of spatial orbitals
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGJacobiOrbitalOptimizer::AP1roGJacobiOrbitalOptimizer(const size_t N_P, const size_t K, const double convergence_threshold, const size_t maximum_number_of_iterations) :
    AP1roGJacobiOrbitalOptimizer(AP1roGGeminalCoefficients(N_P, K), convergence_threshold, maximum_number_of_iterations)
{}


/**
 *  @param G                                the initial geminal coefficients
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
 */
AP1roGJacobiOrbitalOptimizer::AP1roGJacobiOrbitalOptimizer(const AP1roGGeminalCoefficients& G, const double convergence_threshold, const size_t maximum_number_of_iterations) :
    N_P (G.get_N_P()),
    G (G),
    JacobiOrbitalOptimizer(G.get_K(), convergence_threshold, maximum_number_of_iterations)
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
 * 
 *  In the case of this uncoupled AP1roG Jacobi orbital optimizer, we should solve the AP1roG PSEs at the start at every iteration, using the current orbitals
 */
void AP1roGJacobiOrbitalOptimizer::prepareJacobiSpecificConvergenceChecking(const SQHamiltonian<double>& ham_par) {
    
    AP1roGPSESolver pse_solver (this->N_P, ham_par, this->G);
    pse_solver.solve();
    this->G = pse_solver.get_geminal_coefficients();
    this->E = pse_solver.get_electronic_energy();
}


/**
 *  Calculate the trigoniometric polynomial coefficients for the given Jacobi rotation indices
 *
 *  @param p            the index of spatial orbital 1
 *  @param q            the index of spatial orbital 2
 */
void AP1roGJacobiOrbitalOptimizer::calculateJacobiCoefficients(const SQHamiltonian<double>& ham_par, const size_t p, const size_t q) {

    const size_t K = this->dim;
    const auto& h = ham_par.core().parameters();
    const auto& g = ham_par.twoElectron().parameters();
    const auto& G = this->G;

    // Implementation of the Jacobi rotation coefficients with disjoint cases for p and q

    // Occupied-occupied rotations: if p <= N_P and q <= N_P for computers
    if ((p < this->N_P) && (q < this->N_P)) {

        this->A1 = 0.0;
        this->B1 = 0.0;
        this->C1 = 0.0;


        for (size_t b = this->N_P; b < K; b++) {
            this->A1 -= 0.5 * (g(b,p,b,p) - g(b,q,b,q)) * (G(p,b) - G(q,b));
            this->B1 += 0.5 * (g(b,p,b,p) - g(b,q,b,q)) * (G(p,b) - G(q,b));
            this->C1 += g(b,p,b,q) * (G(q,b) - G(p,b));
        }
    }


    // Occupied-virtual rotations: if p > N_P and q <= N_P for computers
    else if ((p >= this->N_P) && (q < this->N_P)) {

        this->A2 = 0.0;
        this->B2 = 0.0;
        this->C2 = 0.0;
        this->D2 = 0.0;
        this->E2 = 0.0;


        this->A2 += h(p,p) - h(q,q) + 0.375 * (g(p,p,p,p) + g(q,q,q,q)) * (1 - G(q,p)) - 0.25 * g(p,p,q,q) * (7 + G(q,p)) + 0.5 * g(p,q,p,q) * (3 + G(q,p));
        this->B2 += h(q,q) - h(p,p) + 2 * g(p,p,q,q) + 0.5 * (g(p,p,p,p) + g(q,q,q,q)) * (G(q,p) - 1) - g(p,q,p,q) * (1 + G(q,p));
        this->C2 += 2 * h(p,q) + (g(p,p,p,q) - g(p,q,q,q)) * (1 - G(q,p));
        this->D2 += 0.125 * (g(p,p,p,p) + g(q,q,q,q) - 2 * (g(p,p,q,q) + 2 * g(p,q,p,q))) * (1 - G(q,p));
        this->E2 += 0.5 * (g(p,p,p,q) - g(p,q,q,q)) * (G(q,p) - 1);

        for (size_t j = 0; j < this->N_P; j++) {
            this->A2 += 2 * (g(j,j,p,p) - g(j,j,q,q)) - 0.5 * (g(j,p,j,p) - g(j,q,j,q)) * (2 + G(j,p));
            this->B2 += 2 * (g(j,j,q,q) - g(j,j,p,p)) + 0.5 * (g(j,p,j,p) - g(j,q,j,q)) * (2 + G(j,p));
            this->C2 += 4 * g(j,j,p,q) - g(j,p,j,q) * (2 + G(j,p));
        }

        for (size_t b = this->N_P; b < K; b++) {
            this->A2 += 0.5 * (g(b,p,b,p) - g(b,q,b,q)) * G(q,b);
            this->B2 += 0.5 * (g(b,q,b,q) - g(b,p,b,p)) * G(q,b);
            this->C2 += g(b,p,b,q) * G(q,b);
        }
    }


    // Virtual-virtual rotations: if p > N_P and q > N_P for computers
    else if ((p >= this->N_P) && (q >= this->N_P )) {

        this->A3 = 0.0;
        this->B3 = 0.0;
        this->C3 = 0.0;


        for (size_t j = 0; j < this->N_P; j++) {
            this->A3 -= 0.5 * (g(j,p,j,p) - g(j,q,j,q)) * (G(j,p) - G(j,q));
            this->B3 += 0.5 * (g(j,p,j,p) - g(j,q,j,q)) * (G(j,p) - G(j,q));
            this->C3 += g(j,p,j,q) * (G(j,q) - G(j,p));
        }
    }

    else {  // this means that p <= q
        throw std::invalid_argument("AP1roGJacobiOrbitalOptimizer::calculateJacobiCoefficients(const SQHamiltonian<double>&, const size_t, const size_t): The given p and q are invalid: p must be larger than q.");
    }
}


/**
 *  @param ham_par      the current Hamiltonian parameters
 *  @param p            the index of spatial orbital 1
 *  @param q            the index of spatial orbital 2
 *
 *  @return the angle for which the derivative of the scalar function after the Jacobi rotation is zero (and the second derivative is positive), using the current trigoniometric polynomial coefficients
 */
double AP1roGJacobiOrbitalOptimizer::calculateOptimalRotationAngle(const SQHamiltonian<double>& ham_par, const size_t p, const size_t q) const {

    // Implementation of the optimal Jacobi rotation angle with disjoint cases for p and q

    // Occupied-occupied rotations: if p <= N_P and q <= N_P for computers
    if ((p < this->N_P) && (q < this->N_P)) {
        const double denominator = std::sqrt(std::pow(this->B1, 2) + std::pow(this->C1, 2));

        // If the denominator is almost zero, the Jacobi rotation is redundant: the corresponding angle of a 'non'-rotation is 0.0
        if (denominator < 1.0e-08) {
            return 0.0;
        }
        return 0.5 * std::atan2(-this->C1 / denominator, -this->B1 / denominator);  // std::atan2(y,x) = tan^-1(y/x)
    }

    // Occupied-virtual rotations: if p > N_P and q <= N_P for computers
    else if ((p >= this->N_P) && (q < this->N_P)) {

        const auto& cmp = this->comparer();
        std::priority_queue<pair_type, std::vector<pair_type>, decltype(cmp)> queue (cmp);

        // Construct a lambda gradient function
        VectorFunction gradient_function = [this] (const VectorX<double>& x) {
            VectorX<double> gradient_vec (1);
            gradient_vec << (-2*this->B2 * std::sin(2*x(0)) + 2*this->C2 * std::cos(2*x(0)) - 4*this->D2 * std::sin(4*x(0)) + 4*this->E2 * std::cos(4*x(0)));
            return gradient_vec;
        };

        // Construct a lambda Hessian function
        MatrixFunction hessian_function = [this] (const VectorX<double>& x) {
            SquareMatrix<double> hessian_matrix (1);
            hessian_matrix << (-4*this->B2 * std::cos(2*x(0)) - 2*this->C2 * std::sin(2*x(0)) - 16*this->D2 * std::cos(4*x(0)) - 16*this->E2 * std::sin(4*x(0)));
            return hessian_matrix;
        };


        // Use three initial guesses to get the minimum
        const auto half_pi = boost::math::constants::half_pi<double>();
        const double quarter_pi = half_pi / 2;
        const std::vector<double> theta_values {0.0, half_pi, quarter_pi};
        for (const auto& theta : theta_values) {
            VectorX<double> theta_vec (1);  // we can't implicitly convert a float to an VectorX<double> so we make it ourselves
            theta_vec << theta;

            NewtonMinimizer minimizer (theta_vec, gradient_function, hessian_function);
            minimizer.solve();
            const double theta_min = minimizer.get_solution()(0);  // get inside the VectorX<double>
            JacobiRotationParameters jacobi_rot_par {p, q, theta_min};

            const double E_change = this->calculateScalarFunctionChange(ham_par, jacobi_rot_par);

            queue.emplace(jacobi_rot_par, E_change);  // construct a pair_type
        }  // for theta

        const double optimal_theta = queue.top().first.get_angle();

        VectorX<double> theta_min_vec (1);  // we can't implicitly convert a float to an VectorX<double> so we make it ourselves
        theta_min_vec << optimal_theta;;

        assert(hessian_function(theta_min_vec)(0,0) > 0);  // the Hessian of the minimal value of the three must be positive, otherwise we're not in a minimum

        return optimal_theta;
    }


    // Virtual-virtual rotations: if p > N_P and q > N_P for computers
    else if ((p >= this->N_P) && (q >= this->N_P )) {
        const double denominator = std::sqrt(std::pow(this->B3, 2) + std::pow(this->C3, 2));

        // If the denominator is almost zero, the Jacobi rotation is redundant: the corresponding angle of a 'non'-rotation is 0.0
        if (denominator < 1.0e-08) {
            return 0.0;
        }
        return 0.5 * std::atan2(-this->C3 / denominator, -this->B3 / denominator);  // std::atan2(y,x) = tan^-1(y/x)
    }


    else {  // this means that p <= q
        throw std::invalid_argument("AP1roGJacobiOrbitalOptimizer::calculateJacobiCoefficients(const SQHamiltonian<double>&, const size_t, const size_t): The given p and q are invalid: p must be larger than q.");
    }
}


/**
 *  @param ham_par              the current Hamiltonian parameters
 *  @param jacobi_rot_par       the Jacobi rotation parameters
 * 
 *  @return the value of the scalar function (i.e. the AP1roG energy) if the given Jacobi rotation parameters would be used to rotate the given Hamiltonian parameters
 */
double AP1roGJacobiOrbitalOptimizer::calculateScalarFunctionChange(const SQHamiltonian<double>& ham_par, const JacobiRotationParameters& jacobi_rot_par) const {

    const size_t p = jacobi_rot_par.get_p();
    const size_t q = jacobi_rot_par.get_q();
    const double theta = jacobi_rot_par.get_angle();


    // I've written everything in terms of cos(2 theta), sin(2 theta), cos(4 theta) and sin(4 theta)
    const double c2 = std::cos(2 * theta);
    const double s2 = std::sin(2 * theta);

    const double c4 = std::cos(4 * theta);
    const double s4 = std::sin(4 * theta);


    // Implementation of the Jacobi rotated energy with disjoint cases for p and q

    // Occupied-occupied rotations: if p <= N_P and q <= N_P for computers
    if ((p < this->N_P) && (q < this->N_P)) {
        return this->A1 + this->B1 * c2 + this->C1 * s2;
    }

    // Occupied-virtual rotations: if p > N_P and q <= N_P for computers
    else if ((p >= this->N_P) && (q < this->N_P)) {
        return this->A2 + this->B2 * c2 + this->C2 * s2 + this->D2 * c4 + this->E2 * s4;
    }

    // Virtual-virtual rotations: if p > N_P and q > N_P for computers
    else if ((p >= this->N_P) && (q >= this->N_P)) {

        return this->A3 + this->B3 * c2 + this->C3 * s2;
    }

    else {  // this means that p <= q
        throw std::invalid_argument("AP1roGJacobiOrbitalOptimizer::calculateJacobiCoefficients(const SQHamiltonian<double>&, const size_t, const size_t): The given p and q are invalid: p must be larger than q.");
    }
}


}  // namespace GQCP
