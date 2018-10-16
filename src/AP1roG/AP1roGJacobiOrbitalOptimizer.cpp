#include "AP1roG/AP1roGJacobiOrbitalOptimizer.hpp"

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
 *  Optimize the AP1roG energy by consequently
 *      - solving the AP1roG equations
 *      - finding the optimal Jacobi transformation (i.e. the one that yields the lowest energy)
 */
void AP1roGJacobiOrbitalOptimizer::orbitalOptimize() {


    GQCG::AP1roGPSESolver pse_solver (2*this->N_P);


    // Solve the PSEs before starting
    this->solvePSE();
    double E_old = this->calculateEnergy();


    bool converged = false;
    constexpr size_t maximum_number_of_iterations = 128;
    size_t iterations = 0;
    while ((!converged) && (iterations < maximum_number_of_iterations)) {

        // Find the Jacobi parameters (p,q,theta) that minimize the energy
        std::priority_queue<JacobiParameters> min_q;  // an ascending queue (on energy) because we have implemented the 'reverse' JacobiParameters::operator<
        min_q.emplace(JacobiParameters {0, 0, 0.0, E_old});  // by including the old energy and invalid Jacobi parameters, we can inform the user if the algorithm didn't find any Jacobi parameters that lower the energy

        for (size_t q = 0; q < this->K; q++) {
            for (size_t p = q+1; p < this->K; p++) {  // loop over p>q
                this->calculateJacobiCoefficients(p,q);

                double theta = this->findOptimalRotationAngle(p, q);
                double E_rotated = this->calculateEnergyAfterJacobiRotation(p, q, theta);

                min_q.emplace(JacobiParameters {p, q, theta, E_rotated});
            }
        }  // loop over p>q

        // Soft check for if there were no Jacobi parameters that lower the energy
        JacobiParameters optimal_jacobi_parameters = min_q.top();
        if (optimal_jacobi_parameters.p == 0 && optimal_jacobi_parameters.q == 0) {
            std::cerr << "We didn't find any Jacobi rotation that lowers the energy." << std::endl;
        }


        // Using the found Jacobi parameters, rotate the basis with the corresponding orthogonal Jacobi matrix
        this->so_basis.rotateJacobi(optimal_jacobi_parameters.p, optimal_jacobi_parameters.q, optimal_jacobi_parameters.theta);


        // Solve the PSEs in the rotated spatial orbital basis
        this->solvePSE();


        // Check for convergence
        double E_new = this->calculateEnergy();
        if (std::abs(E_new - E_old) < threshold) {
            converged = true;
        } else {
            iterations++;
            E_old = E_new;

            if (iterations == maximum_number_of_iterations) {
                throw std::runtime_error("AP1roG::orbitalOptimize: The orbital optimization procedure did not converge.");
            }
        }
    }  // while not converged
}





}  // namespace GQCG
