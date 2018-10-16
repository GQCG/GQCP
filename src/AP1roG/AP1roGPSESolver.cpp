#include "AP1roG/AP1roGPSESolver.hpp"
#include <numopt.hpp>


namespace GQCG {


/*
 * CONSTRUCTORS
 */
/**
 *  Constructor based on a given number of electrons @param N and Hamiltonian parameters @param ham_par
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGPSESolver::AP1roGPSESolver(size_t N, const GQCG::HamiltonianParameters& ham_par) :
    K (ham_par.K),
    ham_par (ham_par),
    N_P (N / 2),
    initial_geminal_coefficients (GQCG::AP1roGGeminalCoefficients(this->N_P, this->K))
{
    // Check if we have an even number of electrons
    if ((N % 2) != 0) {
        throw std::invalid_argument("The given number of electrons is odd.");
    }
}


/**
 *  Constructor based on a given @param molecule and Hamiltonian parameters @param ham_par
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGPSESolver::AP1roGPSESolver(const GQCG::Molecule& molecule, const GQCG::HamiltonianParameters& ham_par) :
    AP1roGPSESolver(molecule.N, ham_par)
{}



/*
 *  PUBLIC METHODS
 */
/**
 *  Calculate the Jacobian element with compound indices (i,a) and (k,c) at the given geminal coefficients @param G
 *
 *      i and k are subscripts, a and c are superscripts
 */
double AP1roGPSESolver::calculateJacobianElement(const AP1roGGeminalCoefficients& G, size_t i, size_t a, size_t k, size_t c) const {

    Eigen::MatrixXd h_SO = this->ham_par.h.get_matrix_representation();
    Eigen::Tensor<double, 4> g_SO = this->ham_par.g.get_matrix_representation();

    double j_el = 0.0;


    // KISS implementation of the calculation of Jacobian elements
    if (i != k) {

        if (a != c) {  // i!=k and a!=c
            return 0.0;
        }

        else {  // i!=k and a == c
            j_el += g_SO(k,i,k,i) - g_SO(k,a,k,a) * G(i,a);

            for (size_t b = this->N_P; b < this->K; b++) {
                if (b != a) {
                    j_el += g_SO(k,b,k,b) * G(i,b);
                }
            }

        }
    }

    else {  // i==k

        if (a != c) {  // i==k and a!=c
            j_el += g_SO(a,c,a,c) - g_SO(i,c,i,c) * G(i,a);

            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    j_el += g_SO(j,c,j,c) * G(j,a);
                }
            }
        }

        else {  // i==k and a==c
            j_el += -2 * g_SO(a,i,a,i) * G(i,a);

            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    j_el += 2 * (2 * g_SO(a,a,j,j) - g_SO(a,j,j,a)) - (2 * g_SO(i,i,j,j) - g_SO(i,j,j,i));
                }
            }

            j_el += 2 * (h_SO(a,a) - h_SO(i,i));

            j_el += g_SO(a,a,a,a) - g_SO(i,i,i,i);

            for (size_t b = this->N_P; b < this->K; b++) {
                if (b != a) {
                    j_el += - g_SO(i,b,i,b) * G(i,b);
                }
            }

            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    j_el += - g_SO(j,a,j,a) * G(j,a);
                }
            }
        }

    }

    return j_el;
}


/**
 *  Calculate and return the Jacobian at the given geminal coefficients @param g
 */
Eigen::MatrixXd AP1roGPSESolver::calculateJacobian(const Eigen::VectorXd& g) const {

    GQCG::AP1roGGeminalCoefficients G (g, this->N_P, this->K);
    size_t number_of_geminal_coefficients = AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K);

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(number_of_geminal_coefficients, number_of_geminal_coefficients);
    // Loop over all Jacobian elements to construct it
    for (size_t mu = 0; mu < number_of_geminal_coefficients; mu++) {
        for (size_t nu = 0; nu < number_of_geminal_coefficients; nu++) {

            // Convert the vector indices mu and nu into matrix indices
            size_t i = G.matrixIndexMajor(nu);
            size_t a = G.matrixIndexMinor(nu);
            size_t k = G.matrixIndexMajor(mu);
            size_t c = G.matrixIndexMinor(mu);

            J(mu, nu) = this->calculateJacobianElement(G, i, a, k, c);
        }
    }

    return J;
}


/**
 *  Calculate the coordinate function at the given geminal coefficients @param G, with given indices.
 *
 *      i is the subscript and a is the superscript
 */
double AP1roGPSESolver::calculateCoordinateFunction(const GQCG::AP1roGGeminalCoefficients& G, size_t i, size_t a) const {

    Eigen::MatrixXd h_SO = this->ham_par.h.get_matrix_representation();
    Eigen::Tensor<double, 4> g_SO = this->ham_par.g.get_matrix_representation();

    double f = 0.0;

    // A KISS implementation of the AP1roG pSE equations
    f += g_SO(a,i,a,i) * (1 - std::pow(G(i,a), 2));

    for (size_t j = 0; j < this->N_P; j++) {
        if (j != i) {
            f += 2 * ((2 * g_SO(a,a,j,j) - g_SO(a,j,j,a)) - (2 * g_SO(i,i,j,j) - g_SO(i,j,j,i))) * G(i,a);
        }
    }

    f += 2 * (h_SO(a,a) - h_SO(i,i)) * G(i,a);

    f += (g_SO(a,a,a,a) - g_SO(i,i,i,i)) * G(i,a);

    for (size_t b = this->N_P; b < this->K; b++) {
        if (b != a) {
            f += (g_SO(a,b,a,b) - g_SO(i,b,i,b) * G(i,a)) * G(i,b);
        }
    }

    for (size_t j = 0; j < this->N_P; j++) {
        if (j != i) {
            f += (g_SO(j,i,j,i) - g_SO(j,a,j,a) * G(i,a)) * G(j,a);
        }
    }

    for (size_t b = this->N_P; b < this->K; b++) {
        if (b != a) {

            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    f += g_SO(j,b,j,b) * G(j,a) * G(i,b);
                }
            }

        }
    }

    return f;
}


/**
 *  Calculate the coordinate functions for the pSEs at the given geminal coefficients @param g. This returns a vector F in which every entry is one of the coordinate functions
 */
Eigen::VectorXd AP1roGPSESolver::calculateCoordinateFunctions(const Eigen::VectorXd& g) const {

    GQCG::AP1roGGeminalCoefficients G (g, this->N_P, this->K);
    size_t number_of_geminal_coefficients = AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K);

    Eigen::VectorXd F = Eigen::VectorXd::Zero(number_of_geminal_coefficients);  // the vector of coordinate functions
    // Loop over all the F elements to construct it
    for (size_t mu = 0; mu < number_of_geminal_coefficients; mu++) {

        // Convert the vector indices mu into matrix indices
        size_t i = G.matrixIndexMajor(mu);
        size_t a = G.matrixIndexMinor(mu);

        F(mu) = this->calculateCoordinateFunction(G, i, a);
    }

    return F;
}


/**
 *  Set up and solve the projected Schrödinger equations for AP1roG
 */
void AP1roGPSESolver::solve() {

    // Solve the AP1roG equations using a Newton-based algorithm

    numopt::VectorFunction f = [this](const Eigen::VectorXd& x) { return this->calculateCoordinateFunctions(x); };
    numopt::JacobianFunction J = [this](const Eigen::VectorXd& x) { return this->calculateJacobian(x); };


    Eigen::VectorXd x0 = this->initial_geminal_coefficients.asVector();
    numopt::syseq::NewtonSystemOfEquationsSolver syseq_solver (x0, f, J);
    syseq_solver.solve();


    // Set the solution
    GQCG::AP1roGGeminalCoefficients geminal_coefficients (syseq_solver.get_solution(), this->N_P, this->K);
    double electronic_energy = GQCG::calculateAP1roGEnergy(geminal_coefficients, this->ham_par);
    this->solution = GQCG::AP1roG(geminal_coefficients, electronic_energy);
}


}  // namespace GQCG
