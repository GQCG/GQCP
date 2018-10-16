#include "AP1roG/AP1roGJacobiOrbitalOptimizer.hpp"

#include <cmath>
#include <queue>

#include <boost/math/constants/constants.hpp>
#include <numopt.hpp>

#include "AP1roG/AP1roGPSESolver.hpp"


namespace GQCG {


/*
 * CONSTRUCTORS
 */
/**
 *  Constructor based on a given number of electron pairs @param N_P, Hamiltonian parameters @param ham_par, a threshold for the orbital optimization @param oo_threshold and a @param maximum_number_of_oo_iterations
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGJacobiOrbitalOptimizer::AP1roGJacobiOrbitalOptimizer(size_t N_P, const GQCG::HamiltonianParameters& ham_par, double oo_threshold, const size_t maximum_number_of_oo_iterations) :
    K (ham_par.K),
    ham_par (ham_par),
    N_P (N_P),
    oo_threshold (oo_threshold),
    maximum_number_of_oo_iterations (maximum_number_of_oo_iterations)
{}


/**
 *  Constructor based on a given @param molecule, Hamiltonian parameters @param ham_par, a threshold for the orbital optimization @param oo_threshold and a @param maximum_number_of_oo_iterations
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGJacobiOrbitalOptimizer::AP1roGJacobiOrbitalOptimizer(const GQCG::Molecule& molecule, const GQCG::HamiltonianParameters& ham_par, double oo_threshold, const size_t maximum_number_of_oo_iterations) :
    AP1roGJacobiOrbitalOptimizer(molecule.N/2, ham_par, oo_threshold, maximum_number_of_oo_iterations)
{
    // Check if we have an even number of electrons
    if ((molecule.N % 2) != 0) {
        throw std::invalid_argument("The given number of electrons is odd.");
    }
}



/*
 *  PUBLIC METHODS
 */
/**
 *  Given the two indices of spatial orbitals @param p and @param q that will be Jacobi-rotated, and the geminal coefficients @param G,     calculate the coefficients (which are @members)
 *      - A1, B1, C1            to be used in occupied-occupied rotations
 *      - A2, B2, C2, D2, E2    to be used in occupied-virtual rotations
 *      - A3, B3, C3            to be used in virtual-virtual rotations
 */
void AP1roGJacobiOrbitalOptimizer::calculateJacobiCoefficients(size_t p, size_t q, const GQCG::AP1roGGeminalCoefficients& G) {

    Eigen::MatrixXd h_SO = this->ham_par.h.get_matrix_representation();
    Eigen::Tensor<double, 4> g_SO = this->ham_par.g.get_matrix_representation();


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
 *  Calculate the AP1roG energy given the geminal coefficients @param G after the application of a Jacobi rotation with the parameters @param jacobi_rotation_parameters
 */
double AP1roGJacobiOrbitalOptimizer::calculateEnergyAfterJacobiRotation(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters, const GQCG::AP1roGGeminalCoefficients& G) const {

    if (!this->are_calculated_jacobi_coefficients) {
        throw std::runtime_error("calculateEnergyAfterJacobiRotation: You haven't calculated the Jacobi coefficients yet. You should call AP1roG::calculateJacobiCoefficients before calling this function.");
    }
    // Note that this throw only warns the first time: it does not detect if there have been changes to the coefficients, for example after a recalculation


    size_t p = jacobi_rotation_parameters.get_p();
    size_t q = jacobi_rotation_parameters.get_q();
    double theta = jacobi_rotation_parameters.get_angle();


    // The formula I have derived is an energy CORRECTION due to the Jacobi rotation, so we initialize the rotated energy by the initial energy
    double E = GQCG::calculateAP1roGEnergy(G, this->ham_par);

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
 *  Given a Jacobi pair @param p and @param q and the geminal coefficients @param G, @return the optimal rotation angle, i.e. the angle for which the derivative
 *  of the energy after the Jacobi rotation is zero (and the second derivative is positive).
 */
double AP1roGJacobiOrbitalOptimizer::findOptimalRotationAngle(size_t p, size_t q, const GQCG::AP1roGGeminalCoefficients& G) const {

    if (!this->are_calculated_jacobi_coefficients) {
        throw std::runtime_error("findOptimalRotationAngle: You haven't calculated the Jacobi coefficients yet. You should call AP1roG::calculateJacobiCoefficients before calling this function.");
    }
    // Note that this throw only warns the first time: it does not detect if there have been changes to the coefficients, for example after a recalculation


    // Implementation of the optimal Jacobi rotation angle with disjoint cases for p and q

    // Occupied-occupied rotations: if p <= N_P and q <= N_P for computers
    if ((p < this->N_P) && (q < this->N_P)) {
        double denominator = std::sqrt(std::pow(this->B1, 2) + std::pow(this->C1, 2));
        return 0.5 * std::atan2(-this->B1 / denominator, -this->C1 / denominator);
    }


    // Occupied-virtual rotations: if p > N_P and q <= N_P for computers
    else if ((p >= this->N_P) && (q < this->N_P)) {

        std::priority_queue<JacobiRotationEnergy> min_q;  // an ascending queue (on energy) because we have implemented the 'reverse' JacobiParameters::operator<

        // Construct a lambda gradient function
        numopt::GradientFunction gradient_function = [this](const Eigen::VectorXd& x) {
            Eigen::VectorXd gradient_vec (1);
            gradient_vec << (-2*this->B2 * std::sin(2*x(0)) + 2*this->C2 * std::cos(2*x(0)) - 4*this->D2 * std::sin(4*x(0)) + 4*this->E2 * std::cos(4*x(0)));
            return gradient_vec;
        };

        // Construct a lambda Hessian function
        numopt::HessianFunction hessian_function = [this](const Eigen::VectorXd& x) {
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

            numopt::minimization::NewtonMinimizer minimizer (theta_vec, gradient_function, hessian_function);
            minimizer.solve();
            double theta_min = minimizer.get_solution()(0);  // get inside the Eigen::VectorXd

            GQCG::JacobiRotationParameters jacobi_rot_par {p, q, theta_min};

            double E_min = this->calculateEnergyAfterJacobiRotation(jacobi_rot_par, G);
            min_q.emplace(JacobiRotationEnergy {jacobi_rot_par, E_min});
        }  // for theta

        Eigen::VectorXd theta_min_vec (1);  // we can't implicitly convert a float to an Eigen::VectorXd so we make it ourselves
        theta_min_vec << min_q.top().jacobi_rotation_parameters.angle;

        assert(hessian_function(theta_min_vec)(0,0) > 0);  // the Hessian of the minimal value of the three must be positive, otherwise we're not in a minimum

        return min_q.top().jacobi_rotation_parameters.angle;
    }


    // Virtual-virtual rotations: if p > N_P and q > N_P for computers
    else if ((p >= this->N_P) && (q >= this->N_P )) {
        double denominator = std::sqrt(std::pow(this->B3, 2) + std::pow(this->C3, 2));
        return 0.5 * std::atan2(-this->B3 / denominator, -this->C3 / denominator);
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
    GQCG::AP1roGPSESolver pse_solver (this->N_P, this->ham_par);
    pse_solver.solve();
    auto G = pse_solver.get_solution().get_geminal_coefficients();
    std::cout << "G as vector: " << std::endl << G.asVector() << std::endl << std::endl;
    double E_old = GQCG::calculateAP1roGEnergy(G, this->ham_par);
    std::cout << "E_old: " << E_old << std::endl;


    size_t iterations = 0;
    while (!(this->is_converged)) {

        // Find the Jacobi parameters (p,q,theta) that minimize the energy
        std::priority_queue<JacobiRotationEnergy> min_q;  // an ascending queue (on energy) because we have implemented the 'reverse' JacobiRotationEnergy::operator<

        for (size_t q = 0; q < this->K; q++) {
            for (size_t p = q+1; p < this->K; p++) {  // loop over p>q
                this->calculateJacobiCoefficients(p, q, G);

                double theta = this->findOptimalRotationAngle(p, q, G);
                GQCG::JacobiRotationParameters jacobi_rot_par {p, q, theta};
                double E_rotated = this->calculateEnergyAfterJacobiRotation(jacobi_rot_par, G);

                min_q.emplace(JacobiRotationEnergy {jacobi_rot_par, E_rotated});
            }
        }  // loop over p>q

        auto optimal_jacobi_parameters = min_q.top().jacobi_rotation_parameters;

        // Check for if there were no Jacobi parameters that lower the energy
        // TODO


        // Using the found Jacobi parameters, rotate the basis with the corresponding orthogonal Jacobi matrix
        this->ham_par.rotate(optimal_jacobi_parameters);


        // Solve the PSEs in the rotated spatial orbital basis
        // UPDATE INITIAL GUESS
        GQCG::AP1roGPSESolver pse_solver (this->N_P, this->ham_par);
        pse_solver.solve();
        auto G = pse_solver.get_solution().get_geminal_coefficients();
        double E_new = GQCG::calculateAP1roGEnergy(G, this->ham_par);
        std::cout << "E_new: " << E_new << std::endl;

        // Check for convergence
        if (std::abs(E_new - E_old) < this->oo_threshold) {
            this->is_converged = true;

            // Set the solution
            double electronic_energy = GQCG::calculateAP1roGEnergy(G, this->ham_par);
            this->solution = GQCG::AP1roG(G, electronic_energy);
        } else {
            iterations++;
            E_old = E_new;

            if (iterations == this->maximum_number_of_oo_iterations) {
                throw std::runtime_error("The orbital optimization procedure did not converge.");
            }
        }
    }  // while not converged
}


}  // namespace GQCG
