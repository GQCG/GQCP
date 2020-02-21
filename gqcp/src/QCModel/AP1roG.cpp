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
#include "QCModel/Geminals/AP1roG.hpp"


namespace GQCP {


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param G                    the converged AP1roG geminal coefficients
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the AP1roG electronic energy
 */
double QCModel::AP1roG::calculateEnergy(const AP1roGGeminalCoefficients& G, const SQHamiltonian<double>& sq_hamiltonian) {

    // Prepare some variables
    const auto& h = sq_hamiltonian.core().parameters();
    const auto& g = sq_hamiltonian.twoElectron().parameters();


    // KISS implementation of the AP1roG energy
    double E = 0.0;
    for (size_t j = 0; j < G.get_N_P(); j++) {
        E += 2 * h(j,j);

        for (size_t k = 0; k < G.get_N_P(); k++) {
            E += 2 * g(k,k,j,j) - g(k,j,j,k);
        }

        for (size_t b = G.get_N_P(); b < G.get_K(); b++) {
            E += g(j,b,j,b) * G(j,b);
        }
    }

    return E;
}


/**
 *  @param G                the AP1roG geminal coefficients
 *  @param multipliers      the AP1roG Lagrangian multipliers
 *
 *  @return the AP1roG response 1-DM
 */
OneRDM<double> QCModel::AP1roG::calculate1RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers) {

    OneRDM<double> D = OneRDM<double>::Zero(G.get_K(), G.get_K());

    size_t N_P = G.get_N_P();
    size_t K = G.get_K();


    // KISS-implementation of the formulas

    // Occupied part
    for (size_t i = 0; i < N_P; i++) {
        double sum {0.0};

        for (size_t a = N_P; a < K; a++) {
            sum += multipliers(i,a) * G(i,a);
        }

        D(i,i) = 2 * (1 - sum);
    }


    // Virtual part
    for (size_t a = N_P; a < K; a++) {
        double sum {0.0};

        for (size_t i = 0; i < N_P; i++) {
            sum += multipliers(i,a) * G(i,a);
        }

        D(a,a) = 2 * sum;
    }

    return D;
}


/**
 *  @param G                the AP1roG geminal coefficients
 *  @param multipliers      the AP1roG Lagrangian multipliers
 *
 *  @return the AP1roG response number 2-RDM (the Delta-matrix in the notes)
 */
SquareMatrix<double> QCModel::AP1roG::calculateNumber2RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers) {

    size_t N_P = G.get_N_P();
    size_t K = G.get_K();

    SquareMatrix<double> Delta = SquareMatrix<double>::Zero(K, K);


    // KISS-implementation
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {

            if ((p < N_P) && (q < N_P)) {  // occupied-occupied block
                size_t i = p;
                size_t j = q;

                double sum {0.0};

                for (size_t a = N_P; a < K; a++) {
                    sum += multipliers(i,a) * G(i,a) + multipliers(j,a) * G(j,a);

                    if (i == j) {
                        sum -= multipliers(i,a) * G(i,a);
                    }
                }

                Delta(i,j) = 4 * (1 - sum);
            }  // occupied-occupied block


            else if ((p >= N_P) && (q >= N_P)) {  // virtual-virtual block
                size_t a = p;
                size_t b = q;

                if (a == b) {
                    double sum {0.0};

                    for (size_t i = 0; i < N_P; i++) {
                        sum += multipliers(i,a) * G(i,a);
                    }

                    Delta(a,b) = 4 * sum;
                }
            }  // virtual-virtual


            else {  // occupied-virtual and virtual-occupied block

                if (p < q) {  // and afterwards set Delta(i,a) = Delta(a,i)
                    size_t i = p;
                    size_t a = q;
                    double sum {0.0};

                    for (size_t j = 0; j < N_P; j++) {
                        if (j != i) {
                            sum += multipliers(j,a) * G(j,a);
                        }
                    }

                    Delta(i,a) = 4 * sum;
                    Delta(a,i) = Delta(i,a);
                }
            }  // occupied-virtual and virtual-occupied block

        }
    }

    return Delta;
}


/**
 *  @param G                the AP1roG geminal coefficients
 *  @param multipliers      the AP1roG Lagrangian multipliers
 *
 *  @return the AP1roG response pair 2-RDM (the Pi-matrix in the notes)
 */
SquareMatrix<double> QCModel::AP1roG::calculatePair2RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers) {

    size_t N_P = G.get_N_P();
    size_t K = G.get_K();

    SquareMatrix<double> Pi = SquareMatrix<double>::Zero(K, K);


    // KISS-implementation
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {

            if ((p < N_P) && (q < N_P)) {  // occupied-occupied block
                size_t i = p;
                size_t j = q;

                double sum {0.0};
                for (size_t a = N_P; a < K; a++) {
                    sum += multipliers(j,a) * G(i,a);
                }


                if (i == j) {  // diagonal occupied part
                    Pi(i,j) += 1.0 - sum;
                } else {
                    Pi(i,j) += sum;
                }

            }  // occupied-occupied block


            else if ((p < N_P) && (q >= N_P)) {  // occupied-virtual block
                size_t i = p;
                size_t a = q;

                double first_sum {0.0};
                for (size_t j = 0; j < N_P; j++) {
                    for (size_t b = N_P; b < K; b++) {
                        first_sum += multipliers(j,b) * G(j,b);
                    }
                }


                double second_sum {0.0};
                for (size_t j = 0; j < N_P; j++) {
                    for (size_t b = N_P; b < K; b++) {
                        if ((j != i) && (b != a)) {
                            second_sum += multipliers(j,b) * (G(i,a) * G(j,b) + G(j,a) * G(i,b));
                        }
                    }
                }

                Pi(i,a) = (1 - first_sum) * G(i,a) + second_sum;
            }  // occupied-virtual


            else if ((p >= N_P) && (q < N_P)) {  // virtual-occupied block
                size_t a = p;
                size_t i = q;

                Pi(a,i) = multipliers(i,a);
            }


            else {  // virtual-virtual block
                size_t a = p;
                size_t b = q;

                double sum {0.0};
                for (size_t i = 0; i < N_P; i++) {
                    sum += multipliers(i,a) * G(i,b);
                }

                Pi(a,b) = sum;
            }

        }
    }

    return Pi;
}


/**
 *  @param G                the AP1roG geminal coefficients
 *  @param multipliers      the AP1roG Lagrangian multipliers
 *
 *  @return the AP1roG response 2-DM
 */
TwoRDM<double> QCModel::AP1roG::calculate2RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers) {

    size_t K = G.get_K();
    TwoRDM<double> d (K);
    d.setZero();


    auto Delta = QCModel::AP1roG::calculateNumber2RDM(G, multipliers);
    auto Pi = QCModel::AP1roG::calculatePair2RDM(G, multipliers);

    // KISS-implementation
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                for (size_t s = 0; s < K; s++) {

                    if ((p == r) && (q == s)) {
                        d(p,q,r,s) += 2 * Pi(p,q);
                    }

                    if ((p == q) && (r == s) && (p != r)) {
                        d(p,q,r,s) += Delta(p,r);
                    }

                    if ((p == s) && (q == r) && (p != q)) {
                        d(p,q,r,s) -= 0.5 * Delta(p,q);
                    }


                }
            }
        }
    }  // spatial orbital loops

    return d;

}


}  // namespace GQCP
