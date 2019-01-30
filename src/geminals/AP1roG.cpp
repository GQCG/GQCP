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
#include "geminals/AP1roG.hpp"


namespace GQCP {


/**
 *  @param G            the converged AP1roG geminal coefficients
 *  @param ham_par      Hamiltonian parameters in an orthonormal spatial orbital basis
 *
 *  @return the AP1roG electronic energy
 */
double calculateAP1roGEnergy(const AP1roGGeminalCoefficients& G, const HamiltonianParameters& ham_par) {

    OneElectronOperator h_SO = ham_par.get_h();
    TwoElectronOperator g_SO = ham_par.get_g();


    // KISS implementation of the AP1roG energy
    double E = 0.0;
    for (size_t j = 0; j < G.get_N_P(); j++) {
        E += 2 * h_SO(j,j);

        for (size_t k = 0; k < G.get_N_P(); k++) {
            E += 2 * g_SO(k,k,j,j) - g_SO(k,j,j,k);
        }

        for (size_t b = G.get_N_P(); b < G.get_K(); b++) {
            E += g_SO(j,b,j,b) * G(j,b);
        }
    }

    return E;
}


/**
 *  @param G            the AP1roG geminal coefficients
 *  @param Q            the AP1roG bivariational coefficients
 *
 *  @return the overlap between the bivariational coefficients and the geminal coefficients, i.e. <Phi(q)|Psi(p)>
 */
double calculateOverlap(const AP1roGGeminalCoefficients& G, const BivariationalCoefficients& Q) {

    double overlap = Q.q0;
    AP1roGVariables q = Q.q;

    for (size_t i = 0; i < G.get_N_P(); i++) {
        for (size_t a = G.get_N_P(); a < G.get_K(); a++) {
            overlap += q(i,a) * G(i,a);
        }
    }

    return overlap;
}


/**
 *  @param G            the AP1roG geminal coefficients
 *  @param Q            the AP1roG bivariational coefficients
 *
 *  @return the AP1roG 1-DM
 */
OneRDM calculate1RDM(const AP1roGGeminalCoefficients& G, const BivariationalCoefficients& Q) {

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(G.get_K(), G.get_K());
    double overlap = calculateOverlap(G, Q);

    AP1roGVariables q = Q.q;
    size_t N_P = G.get_N_P();
    size_t K = G.get_K();


    // KISS-implementation of the formulas

    // Occupied part
    for (size_t i = 0; i < N_P; i++) {
        double intermediate = Q.q0;

        for (size_t k = 0; k < N_P; k++) {
            for (size_t a = N_P; a < K; a++) {
                if (k != i) {
                    intermediate += q(k,a) * G(k,a);
                }
            }
        }

        D(i,i) = 2 / overlap * intermediate;
    }


    // Virtual part
    for (size_t a = N_P; a < K; a++) {
        double intermediate = 0.0;

        for (size_t i = 0; i < N_P; i++) {
            intermediate += q(i,a) * G(i,a);
        }

        D(a,a) = 2 / overlap * intermediate;
    }

    return OneRDM(D);
}



/**
 *  @param G            the AP1roG geminal coefficients
 *  @param Q            the AP1roG bivariational coefficients
 *
 *  @return the AP1roG number 2-RDM (the Delta-matrix in the notes)
 */
Eigen::MatrixXd calculateNumber2RDM(const AP1roGGeminalCoefficients& G, const BivariationalCoefficients& Q) {

    size_t N_P = G.get_N_P();
    size_t K = G.get_K();
    double overlap = calculateOverlap(G, Q);

    Eigen::MatrixXd Delta = Eigen::MatrixXd::Zero(K, K);


    // KISS-implementation
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {

            if ((p < N_P) && (q < N_P)) {  // occupied-occupied block
                size_t i = p;
                size_t j = q;

                double intermediate = 0.0;

                for (size_t c = N_P; c < K; c++) {
                    intermediate += Q.q(i,c) * G(i,c) + Q.q(j,c) * G(j,c);

                    if (i == j) {
                        intermediate -= Q.q(i,c) * G(i,c);
                    }
                }

                Delta(i,j) = 4 * (1 - intermediate / overlap);
            }  // occupied-occupied block


            else if ((p >= N_P) && (q >= N_P)) {  // virtual-virtual block
                size_t a = p;
                size_t b = q;

                if (a == b) {
                    double intermediate = 0.0;

                    for (size_t k = 0; k < N_P; k++) {
                        intermediate += Q.q(k,a) * G(k,a);
                    }

                    Delta(a,b) = 4 * intermediate / overlap;
                }
            }  // virtual-virtual


            else {  // occupied-virtual and virtual-occupied block

                if (p < q) {  // and afterwards set Delta(i,a) = Delta(a,i)
                    size_t i = p;
                    size_t a = q;
                    double intermediate = 0.0;

                    for (size_t k = 0; k < N_P; k++) {
                        if (k != i) {
                            intermediate += Q.q(k,a) * G(k,a);
                        }
                    }

                    Delta(i,a) = 4 * intermediate / overlap;
                    Delta(a,i) = Delta(i,a);
                }
            }  // occupied-virtual and virtual-occupied block

        }
    }

    return Delta;
}


/**
 *  @param G            the AP1roG geminal coefficients
 *  @param Q            the AP1roG bivariational coefficients
 *
 *  @return the AP1roG pair 2-RDM (the Pi-matrix in the notes)
 */
Eigen::MatrixXd calculatePair2RDM(const AP1roGGeminalCoefficients& G, const BivariationalCoefficients& Q) {

    size_t N_P = G.get_N_P();
    size_t K = G.get_K();
    double overlap = calculateOverlap(G, Q);

    Eigen::MatrixXd Pi = Eigen::MatrixXd::Zero(K, K);


    // KISS-implementation
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {

            if ((p < N_P) && (q < N_P)) {  // occupied-occupied block
                size_t i = p;
                size_t j = q;

                double intermediate = 0.0;
                for (size_t a = N_P; a < K; a++) {
                    intermediate += Q.q(j,a) * G(i,a);

                    if (i == j) {
                        intermediate -= 2 * Q.q(i,a) * G(i,a);
                    }
                }


                if (i == j) {  // delta_ij
                    Pi(i,j) += 1.0;
                }

                Pi(i,j) += intermediate / overlap;
            }  // occupied-occupied block


            else if ((p < N_P) && (q >= N_P)) {  // occupied-virtual block
                size_t i = p;
                size_t a = q;

                double intermediate = Q.q0 * G(i,a);
                for (size_t j = 0; j < N_P; j++) {
                    for (size_t b = N_P; b < K; b++) {
                        if ((j != i) && (b != a)) {
                            intermediate += Q.q(j,b) * (G(i,a) * G(j,b) + G(j,a) * G(i,b));
                        }
                    }
                }

                Pi(i,a) = intermediate / overlap;
            }  // occupied-virtual


            else if ((p >= N_P) && (q < N_P)) {  // virtual-occupied block
                size_t a = p;
                size_t i = q;

                Pi(a,i) = Q.q(i,a) / overlap;
            }


            else {  // virtual-virtual block
                size_t a = p;
                size_t b = q;

                double intermediate = 0.0;
                for (size_t i = 0; i < N_P; i++) {
                    intermediate += Q.q(i,a) * G(i,b);
                }

                Pi(a,b) = intermediate / overlap;
            }

        }
    }

    return Pi;
}


/**
 *  @param G            the AP1roG geminal coefficients
 *  @param Q            the AP1roG bivariational coefficients
 *
 *  @return the AP1roG 2-DM
 */
TwoRDM calculate2RDM(const AP1roGGeminalCoefficients& G, const BivariationalCoefficients& Q) {

    size_t K = G.get_K();
    Eigen::Tensor<double, 4> d (K, K, K, K);
    d.setZero();


    Eigen::MatrixXd Delta = calculateNumber2RDM(G, Q);
    Eigen::MatrixXd Pi = calculatePair2RDM(G, Q);

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

    return TwoRDM(d);
}


}  // namespace GQCP
