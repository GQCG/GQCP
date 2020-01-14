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
#include "QCMethod/Geminals/AP1roGPSEs.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param sq_hamiltonian       the one- and two-electron integrals in an orthonormal orbital basis
 *  @param N_P                  the number of electron pairs
 */
AP1roGPSEs::AP1roGPSEs(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P) :
    N_P (N_P),
    K (sq_hamiltonian.dimension()),
    sq_hamiltonian (sq_hamiltonian)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param G        the AP1roG geminal coefficients
 *  @param i        the subscript for the coordinate function
 *  @param a        the superscript for the coordinate function
 *
 *  @return the value of the coordinate function with given indices (i,a) at the given geminal coefficients
 */
double AP1roGPSEs::calculateCoordinateFunction(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a) const {

    const auto& h = this->sq_hamiltonian.core().parameters();
    const auto& g = this->sq_hamiltonian.twoElectron().parameters();


    double value = 0.0;
    // A KISS implementation of the AP1roG pSE equations
    value += g(a,i,a,i) * (1 - std::pow(G(i,a), 2));

    for (size_t j = 0; j < this->N_P; j++) {
        if (j != i) {
            value += 2 * ((2 * g(a,a,j,j) - g(a,j,j,a)) - (2 * g(i,i,j,j) - g(i,j,j,i))) * G(i,a);
        }
    }

    value += 2 * (h(a,a) - h(i,i)) * G(i,a);

    value += (g(a,a,a,a) - g(i,i,i,i)) * G(i,a);

    for (size_t b = this->N_P; b < this->K; b++) {
        if (b != a) {
            value += (g(a,b,a,b) - g(i,b,i,b) * G(i,a)) * G(i,b);
        }
    }

    for (size_t j = 0; j < this->N_P; j++) {
        if (j != i) {
            value += (g(j,i,j,i) - g(j,a,j,a) * G(i,a)) * G(j,a);
        }
    }

    for (size_t b = this->N_P; b < this->K; b++) {
        if (b != a) {
            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    value += g(j,b,j,b) * G(j,a) * G(i,b);
                }
            }

        }
    }

    return value;
}


/**
 *  @param G            the geminal coefficients
 * 
 *  @return the PSEs, evaluated at the given geminal coefficients
 */
BlockMatrix<double> AP1roGPSEs::calculateCoordinateFunctions(const AP1roGGeminalCoefficients& G) const {

    BlockMatrix<double> F (0, this->N_P, this->N_P, this->K);  // an occupied-virtual matrix

    // Loop over all the elements of F to construct it
    for (size_t i = 0; i < this->N_P; i++) {
        for (size_t a = N_P; a < this->K; a++) {
            F(i,a) = this->calculateCoordinateFunction(G, i, a);
        }
    }

    return F;
}


/**
 *  @return a callable (i.e. with operator()) expression for the coordinate functions: the accepted VectorX<double> argument should contain the geminal coefficients in a column-major representation
 */
VectorFunction AP1roGPSEs::callableCoordinateFunctions() const {

    VectorFunction callable = [this] (const VectorX<double>& x) {
        auto G = AP1roGGeminalCoefficients::FromColumnMajor(x, this->N_P, this->K);
        return this->calculateCoordinateFunctions(G).asVector(); 
    };

    return callable;
}


/**
 *  @param G        the AP1roG geminal coefficients
 *  @param i        the subscript for the coordinate function
 *  @param a        the superscript for the coordinate function
 *  @param j        the subscript for the geminal coefficient
 *  @param b        the superscript for the geminal coefficient
 *
 *  @return the value of the Jacobian element with compound indices (i,a) and (j,b) at the given geminal coefficients
 */
double AP1roGPSEs::calculateJacobianElement(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a, const size_t j, const size_t b) const {

    const auto& h = this->sq_hamiltonian.core().parameters();
    const auto& g = this->sq_hamiltonian.twoElectron().parameters();


    double value = 0.0;
    // KISS implementation of the calculation of Jacobian elements
    if (i == j) {
        value += g(a,b,a,b) - 2 * g(j,b,j,b) * G(j,a);

        for (size_t k = 0; k < this->N_P; k++) {
            value += g(k,b,k,b) * G(k,a);
        }
    }


    if (a == b) {
        value += g(j,i,j,i) - 2 * g(j,b,j,b) * G(i,b);

        for (size_t c = this->N_P; c < this->K; c++) {
            value += g(j,c,j,c) * G(i,c);
        }
    }


    if ((i == j) && (a == b)) {
        value += 2 * (h(a,a) - h(i,i));

        value -= 2 * (2 * g(a,a,i,i) - g(a,i,i,a));

        for (size_t k = 0; k < this->N_P; k++) {
            value += 2 * (2 * g(k,k,a,a) - g(a,k,k,a)) - 2 * (2 * g(i,i,k,k) - g(i,k,k,i));
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

    return value;
}


/**
 *  @param G            the geminal coefficients
 * 
 *  @return the Jacobian, J_{ia,jb} of the PSEs, evaluated at the given geminal coefficients
 */
BlockRankFourTensor<double> AP1roGPSEs::calculateJacobian(const AP1roGGeminalCoefficients& G) const {

    BlockRankFourTensor<double> J (0, this->N_P, this->N_P, this->K, 
                                   0, this->N_P, this->N_P, this->K);  // an occupied-virtual, occupied-virtual tensor


    // Loop over all elements of J to construct it
    for (size_t i = 0; i < this->N_P; i++) {
        for (size_t a = N_P; a < this->K; a++) {
            for (size_t j = 0; j < this->N_P; j++) {
                for (size_t b = N_P; b < this->K; b++) {
                    J(i,a,j,b) = this->calculateJacobianElement(G, i, a, j, b);
                }
            }
        }
    }

    return J;
}


/**
 *  @return a callable (i.e. with operator()) expression for the Jacobian: the accepted VectorX<double> argument should contain the geminal coefficients in a column-major representation
 */
MatrixFunction AP1roGPSEs::callableJacobian() const {

    MatrixFunction callable = [this] (const VectorX<double>& x) { 
        auto G = AP1roGGeminalCoefficients::FromColumnMajor(x, this->N_P, this->K);
        return this->calculateJacobian(G).asMatrix();
    };

    return callable;
}


}  // namespace GQCP
