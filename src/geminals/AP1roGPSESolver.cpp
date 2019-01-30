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

#include "geminals/AP1roG.hpp"
#include "optimization/NewtonSystemOfEquationsSolver.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 */
AP1roGPSESolver::AP1roGPSESolver(size_t N_P, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G) :
    BaseAP1roGSolver(N_P, ham_par, G)
{}


/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGPSESolver::AP1roGPSESolver(size_t N_P, const HamiltonianParameters& ham_par) :
    BaseAP1roGSolver(N_P, ham_par)
{}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 */
AP1roGPSESolver::AP1roGPSESolver(const Molecule& molecule, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G) :
    BaseAP1roGSolver(molecule, ham_par, G)
{
}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGPSESolver::AP1roGPSESolver(const Molecule& molecule, const HamiltonianParameters& ham_par) :
    BaseAP1roGSolver(molecule, ham_par)
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

    OneElectronOperator h_SO = this->ham_par.get_h();
    TwoElectronOperator g_SO = this->ham_par.get_g();

    double j_el = 0.0;


    // KISS implementation of the calculation of Jacobian elements
    if (i != k) {

        if (a != c) {  // i!=k and a!=c
            return 0.0;
        }

        else {  // i!=k and a == c
            j_el += g_SO(k,i,k,i) - 2 * g_SO(k,a,k,a) * G(i,a);

            for (size_t b = this->N_P; b < this->K; b++) {
                j_el += g_SO(k,b,k,b) * G(i,b);
            }

        }
    }

    else {  // i==k

        if (a != c) {  // i==k and a!=c
            j_el += g_SO(a,c,a,c) - 2 * g_SO(i,c,i,c) * G(i,a);

            for (size_t j = 0; j < this->N_P; j++) {
                j_el += g_SO(j,c,j,c) * G(j,a);
            }
        }

        else {  // i==k and a==c

            j_el += 2 * (h_SO(a,a) - h_SO(i,i));

            j_el += g_SO(a,a,a,a) + g_SO(i,i,i,i);

            j_el -= 2 * (2 * g_SO(a,a,i,i) - g_SO(a,i,i,a));


            for (size_t j = 0; j < this->N_P; j++) {
                j_el += 2 * (2 * g_SO(a,a,j,j) - g_SO(a,j,j,a)) - (2 * g_SO(i,i,j,j) - g_SO(i,j,j,i));
            }

            for (size_t j = 0; j < this->N_P; j++) {
                j_el -= g_SO(j,a,j,a) * G(j,a);
            }

            for (size_t b = this->N_P; b < this->K; b++) {
                j_el -= g_SO(i,b,i,b) * G(i,b);
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

    AP1roGGeminalCoefficients G (g, this->N_P, this->K);
    size_t number_of_geminal_coefficients = AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K);

    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(number_of_geminal_coefficients, number_of_geminal_coefficients);
    // Loop over all Jacobian elements to construct it
    for (size_t row_index = 0; row_index < number_of_geminal_coefficients; row_index++) {
        for (size_t column_index = 0; column_index < number_of_geminal_coefficients; column_index++) {

            // In our definitions, we have:
            //      row indices refer to the coordinate functions
            size_t i = G.matrixIndexMajor(row_index);
            size_t a = G.matrixIndexMinor(row_index);

            //      column indices refer to geminal coefficients
            size_t k = G.matrixIndexMajor(column_index);
            size_t c = G.matrixIndexMinor(column_index);

            J(row_index,column_index) = this->calculateJacobianElement(G, i, a, k, c);
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
double AP1roGPSESolver::calculateCoordinateFunction(const AP1roGGeminalCoefficients& G, size_t i, size_t a) const {

    OneElectronOperator h_SO = this->ham_par.get_h();
    TwoElectronOperator g_SO = this->ham_par.get_g();

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

    AP1roGGeminalCoefficients G (g, this->N_P, this->K);
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

    VectorFunction f = [this](const Eigen::VectorXd& x) { return this->calculateCoordinateFunctions(x); };
    MatrixFunction J = [this](const Eigen::VectorXd& x) { return this->calculateJacobian(x); };


    Eigen::VectorXd x0 = this->geminal_coefficients.asVector();
    NewtonSystemOfEquationsSolver syseq_solver (x0, f, J);
    syseq_solver.solve();


    // Set the solution
    this->geminal_coefficients = AP1roGGeminalCoefficients(syseq_solver.get_solution(), this->N_P, this->K);
    this->electronic_energy = calculateAP1roGEnergy(this->geminal_coefficients, this->ham_par);
}


}  // namespace GQCP
