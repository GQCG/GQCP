// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#include "geminals/AP1roGPSESolver.hpp"
#include <numopt.hpp>


namespace GQCP {


/*
 * CONSTRUCTORS
 */
/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 */
AP1roGPSESolver::AP1roGPSESolver(size_t N_P, const GQCP::HamiltonianParameters& ham_par, const GQCP::AP1roGGeminalCoefficients& G) :
    K (ham_par.get_K()),
    ham_par (ham_par),
    N_P (N_P),
    initial_geminal_coefficients (G)
{}

/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGPSESolver::AP1roGPSESolver(size_t N_P, const GQCP::HamiltonianParameters& ham_par) :
    GQCP::AP1roGPSESolver(N_P, ham_par, GQCP::AP1roGGeminalCoefficients(N_P, ham_par.get_K()))
{}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 */
AP1roGPSESolver::AP1roGPSESolver(const GQCP::Molecule& molecule, const GQCP::HamiltonianParameters& ham_par, const GQCP::AP1roGGeminalCoefficients& G) :
    GQCP::AP1roGPSESolver(molecule.get_N()/2, ham_par, G)
{
    // Check if we have an even number of electrons
    if ((molecule.get_N() % 2) != 0) {
        throw std::invalid_argument("The given number of electrons is odd.");
    }
}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGPSESolver::AP1roGPSESolver(const GQCP::Molecule& molecule, const GQCP::HamiltonianParameters& ham_par) :
    AP1roGPSESolver(molecule, ham_par, GQCP::AP1roGGeminalCoefficients(molecule.get_N()/2, ham_par.get_K()))
{}



/*
 *  PUBLIC METHODS
 */
/**
 *  @param G        the AP1roG geminal coefficients
 *  @param i        the subscript for the coordinate function
 *  @param a        the superscript for the coordinate function
 *  @param k        the subscript for the geminal coefficient
 *  @param c        the superscript for the geminal coefficient
 *
 *  @return the Jacobian element with compound indices (i,a) and (k,c) at the given geminal coefficients
 */
double AP1roGPSESolver::calculateJacobianElement(const AP1roGGeminalCoefficients& G, size_t i, size_t a, size_t k, size_t c) const {

    GQCP::OneElectronOperator h_SO = this->ham_par.get_h();
    GQCP::TwoElectronOperator g_SO = this->ham_par.get_g();

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
 *  @param g        the AP1roG geminal coefficients in row-major vector form
 *
 *  @return the Jacobian at the given geminal coefficients
 */
Eigen::MatrixXd AP1roGPSESolver::calculateJacobian(const Eigen::VectorXd& g) const {

    GQCP::AP1roGGeminalCoefficients G (g, this->N_P, this->K);
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
 *  @param G        the AP1roG geminal coefficients
 *  @param i        the subscript for the coordinate function
 *  @param a        the superscript for the coordinate function
 *
 *  @return the coordinate function with given indices (i,a) at the given geminal coefficients
 */
double AP1roGPSESolver::calculateCoordinateFunction(const GQCP::AP1roGGeminalCoefficients& G, size_t i, size_t a) const {

    GQCP::OneElectronOperator h_SO = this->ham_par.get_h();
    GQCP::TwoElectronOperator g_SO = this->ham_par.get_g();

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
 *  @param g        the AP1roG geminal coefficients in row-major vector form
 *
 *  @return the vector of coordinate functions at the given geminal coefficients
 */
Eigen::VectorXd AP1roGPSESolver::calculateCoordinateFunctions(const Eigen::VectorXd& g) const {

    GQCP::AP1roGGeminalCoefficients G (g, this->N_P, this->K);
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
    numopt::MatrixFunction J = [this](const Eigen::VectorXd& x) { return this->calculateJacobian(x); };


    Eigen::VectorXd x0 = this->initial_geminal_coefficients.asVector();
    numopt::syseq::NewtonSystemOfEquationsSolver syseq_solver (x0, f, J);
    syseq_solver.solve();


    // Set the solution
    GQCP::AP1roGGeminalCoefficients geminal_coefficients (syseq_solver.get_solution(), this->N_P, this->K);
    double electronic_energy = GQCP::calculateAP1roGEnergy(geminal_coefficients, this->ham_par);
    this->solution = GQCP::AP1roG(geminal_coefficients, electronic_energy);
}


}  // namespace GQCP
