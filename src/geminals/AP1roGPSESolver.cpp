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
#include "Geminals/AP1roGPSESolver.hpp"

#include "Geminals/AP1roG.hpp"
#include "math/optimization/NewtonSystemOfEquationsSolver.hpp"


namespace GQCP {


/*
 *  PUBLIC METHODS
 */

/**
 *  @param G        the AP1roG geminal coefficients
 *  @param i        the subscript for the coordinate function
 *  @param a        the superscript for the coordinate function
 *  @param j        the subscript for the geminal coefficient
 *  @param b        the superscript for the geminal coefficient
 *
 *  @return the Jacobian element with compound indices (i,a) and (j,b) at the given geminal coefficients
 */
double AP1roGPSESolver::calculateJacobianElement(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a, const size_t j, const size_t b) const {

    const auto& h = this->ham_par.get_h();
    const auto& g = this->ham_par.get_g();

    double value = 0.0;


    // KISS implementation of the calculation of Jacobian elements
    if (i != j) {

        if (a != b) {  // i!=j and a!=b
            return 0.0;
        }

        else {  // i!=j and a == b
            value += g(j,i,j,i) - 2 * g(j,b,j,b) * G(i,b);

            for (size_t c = this->N_P; c < this->K; c++) {
                value += g(j,c,j,c) * G(i,c);
            }

        }
    }

    else {  // i==j

        if (a != b) {  // i==j and a!=b
            value += g(a,b,a,b) - 2 * g(j,b,j,b) * G(j,a);

            for (size_t k = 0; k < this->N_P; k++) {
                value += g(k,b,k,b) * G(k,a);
            }
        }

        else {  // i==j and a==b

            value += 2 * (h(a,a) - h(i,i));

            value -= 2 * (2 * g(a,a,i,i) - g(a,i,i,a));

            for (size_t k = 0; k < this->N_P; k++) {
                value += 2 * (2 * g(k,k,a,a) - g(a,k,k,a)) - (2 * g(i,i,k,k) - g(i,k,k,i));
            }

            for (size_t k = 0; k < this->N_P; k++) {
                if (k != i) {
                    value -= 2 * g(k,a,k,a) * G(k,a);
                }
            }

            for (size_t c = this->N_P; c < this->K; c++) {
                if (c != a) {
                    value -= 2 * g(i,c,i,c) * G(i,c);
                }
            }
        }

    }

    return value;
}


/**
 *  @param G        the AP1roG geminal coefficients
 *
 *  @return the Jacobian at the given geminal coefficients
 */
SquareMatrix<double> AP1roGPSESolver::calculateJacobian(const AP1roGGeminalCoefficients& G) const {

    size_t number_of_geminal_coefficients = AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K);

    // Loop over all Jacobian elements to construct it
    SquareMatrix<double> J = SquareMatrix<double>::Zero(number_of_geminal_coefficients, number_of_geminal_coefficients);
    for (size_t row_index = 0; row_index < number_of_geminal_coefficients; row_index++) {
        for (size_t column_index = 0; column_index < number_of_geminal_coefficients; column_index++) {

            // Using our convention, the Jacobian is defined as df_j^b/dt_i^a:
            //      Column indices refer to the coordinate functions
            size_t j = G.matrixIndexMajor(row_index);
            size_t b = G.matrixIndexMinor(row_index);

            //      Row indices refer to geminal coefficients
            size_t i = G.matrixIndexMajor(column_index);
            size_t a = G.matrixIndexMinor(column_index);

            J(row_index,column_index) = this->calculateJacobianElement(G, j, b, i, a);
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
double AP1roGPSESolver::calculateCoordinateFunction(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a) const {

    const auto& h = this->ham_par.get_h();
    const auto& g = this->ham_par.get_g();

    double f = 0.0;

    // A KISS implementation of the AP1roG pSE equations
    f += g(a,i,a,i) * (1 - std::pow(G(i,a), 2));

    for (size_t j = 0; j < this->N_P; j++) {
        if (j != i) {
            f += 2 * ((2 * g(a,a,j,j) - g(a,j,j,a)) - (2 * g(i,i,j,j) - g(i,j,j,i))) * G(i,a);
        }
    }

    f += 2 * (h(a,a) - h(i,i)) * G(i,a);

    f += (g(a,a,a,a) - g(i,i,i,i)) * G(i,a);

    for (size_t b = this->N_P; b < this->K; b++) {
        if (b != a) {
            f += (g(a,b,a,b) - g(i,b,i,b) * G(i,a)) * G(i,b);
        }
    }

    for (size_t j = 0; j < this->N_P; j++) {
        if (j != i) {
            f += (g(j,i,j,i) - g(j,a,j,a) * G(i,a)) * G(j,a);
        }
    }

    for (size_t b = this->N_P; b < this->K; b++) {
        if (b != a) {

            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    f += g(j,b,j,b) * G(j,a) * G(i,b);
                }
            }

        }
    }

    return f;
}


/**
 *  @param G        the AP1roG geminal coefficients
 *
 *  @return the vector of coordinate functions at the given geminal coefficients
 */
VectorX<double> AP1roGPSESolver::calculateCoordinateFunctions(const AP1roGGeminalCoefficients& G) const {

    size_t number_of_geminal_coefficients = AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K);

    // Loop over all the elements of F to construct it
    VectorX<double> F = VectorX<double>::Zero(number_of_geminal_coefficients);  // the vector of coordinate functions
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

    VectorFunction f = [this](const VectorX<double>& x) { 
        AP1roGGeminalCoefficients G (x, this->N_P, this->K);
        return this->calculateCoordinateFunctions(G); 
    };
    MatrixFunction J = [this](const VectorX<double>& x) { 
        AP1roGGeminalCoefficients G (x, this->N_P, this->K);
        return this->calculateJacobian(G); 
    };


    VectorX<double> x0 = this->geminal_coefficients.asVector();
    NewtonSystemOfEquationsSolver syseq_solver (x0, f, J, this->convergence_threshold, this->maximum_number_of_iterations);
    syseq_solver.solve();


    // Set the solution
    this->geminal_coefficients = AP1roGGeminalCoefficients(syseq_solver.get_solution(), this->N_P, this->K);
    this->electronic_energy = calculateAP1roGEnergy(this->geminal_coefficients, this->ham_par);
}


}  // namespace GQCP
