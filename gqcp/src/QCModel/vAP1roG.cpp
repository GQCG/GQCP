// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "QCModel/Geminals/vAP1roG.hpp"

#include "QCModel/Geminals/AP1roG.hpp"


namespace GQCP {


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param G                the AP1roG geminal coefficients
 *  @param multipliers      the AP1roG Lagrangian multipliers
 *
 *  @return the AP1roG response 1-DM
 */
OneRDM<double> QCModel::vAP1roG::calculate1RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers) {

    // KISS-implementation of the formulas.
    OneRDM<double> D = OneRDM<double>::Zero(G.get_K(), G.get_K());

    const auto orbital_space = G.orbitalSpace();


    // Occupied part.
    for (const auto& i : orbital_space.occupiedIndices()) {
        double sum {0.0};

        for (const auto& a : orbital_space.virtualIndices()) {
            sum += multipliers(i, a) * G(i, a);
        }

        D(i, i) = 2 * (1 - sum);
    }


    // Virtual part.
    for (const auto& a : orbital_space.virtualIndices()) {
        double sum {0.0};

        for (const auto& i : orbital_space.occupiedIndices()) {
            sum += multipliers(i, a) * G(i, a);
        }

        D(a, a) = 2 * sum;
    }

    return D;
}


/**
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *  @param N_P                  the number of electron pairs
 * 
 *  @return the response force (-F_lambda) that is used to solve the linear equations for the Lagrange multipliers lambda in [k_lambda lambda = -F_lambda]
 */
BlockMatrix<double> QCModel::vAP1roG::calculateMultiplierResponseForce(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P) {

    // Prepare some variables.
    const auto& g = sq_hamiltonian.twoElectron().parameters();

    // Create an occupied-virtual orbital space.
    const auto K = sq_hamiltonian.dimension();  // number of spatial orbitals
    const auto orbital_space = OrbitalSpace::OccupiedVirtual(N_P, K);


    BlockMatrix<double> F_lambda {0, N_P, N_P, K};
    for (const auto& i : orbital_space.occupiedIndices()) {
        for (const auto& a : orbital_space.virtualIndices()) {
            F_lambda(i, a) = -g(i, a, i, a);
        }
    }

    return F_lambda;
}


/**
 *  @param G                    the AP1roG geminal coefficients
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 * 
 *  @return the response force constant (k_lambda) that is used to solve the linear equations for the Lagrange multipliers lambda in [k_lambda lambda = -F_lambda]
 */
MatrixX<double> QCModel::vAP1roG::calculateMultiplierResponseForceConstant(const SQHamiltonian<double>& sq_hamiltonian, const AP1roGGeminalCoefficients& G) {

    const auto N_P = G.numberOfElectronPairs();

    const MatrixX<double> k_lambda = QCModel::AP1roG::calculatePSEJacobian(sq_hamiltonian, G).asMatrix().transpose();
    return k_lambda;
}


/**
 *  @param G                the AP1roG geminal coefficients
 *  @param multipliers      the AP1roG Lagrangian multipliers
 *
 *  @return the AP1roG response number 2-RDM (the Delta-matrix in the notes)
 */
SquareMatrix<double> QCModel::vAP1roG::calculateNumber2RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers) {

    const size_t K = G.get_K();
    const auto orbital_space = G.orbitalSpace();

    SquareMatrix<double> Delta = SquareMatrix<double>::Zero(K, K);


    // KISS-implementation
    for (const auto& p : orbital_space.allIndices()) {
        for (const auto& q : orbital_space.allIndices()) {

            if (orbital_space.isIndexOccupied(p) && orbital_space.isIndexOccupied(q)) {  // occupied-occupied block
                const size_t i = p;
                const size_t j = q;

                double sum {0.0};

                for (const auto& a : orbital_space.virtualIndices()) {
                    sum += multipliers(i, a) * G(i, a) + multipliers(j, a) * G(j, a);

                    if (i == j) {
                        sum -= multipliers(i, a) * G(i, a);
                    }
                }

                Delta(i, j) = 4 * (1 - sum);
            }  // occupied-occupied block


            else if (orbital_space.isIndexVirtual(p) && orbital_space.isIndexVirtual(q)) {  // virtual-virtual block
                const size_t a = p;
                const size_t b = q;

                if (a == b) {
                    double sum {0.0};

                    for (const auto& i : orbital_space.occupiedIndices()) {
                        sum += multipliers(i, a) * G(i, a);
                    }

                    Delta(a, b) = 4 * sum;
                }
            }  // virtual-virtual


            else {  // occupied-virtual and virtual-occupied block

                if (p < q) {  // and afterwards set Delta(i,a) = Delta(a,i)
                    const size_t i = p;
                    const size_t a = q;
                    double sum {0.0};

                    for (const auto& j : orbital_space.occupiedIndices()) {
                        if (j != i) {
                            sum += multipliers(j, a) * G(j, a);
                        }
                    }

                    Delta(i, a) = 4 * sum;
                    Delta(a, i) = Delta(i, a);
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
 *  @return the vAP1roG response pair 2-RDM (the Pi-matrix in the notes)
 */
SquareMatrix<double> QCModel::vAP1roG::calculatePair2RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers) {

    const size_t K = G.get_K();
    const auto orbital_space = G.orbitalSpace();


    SquareMatrix<double> Pi = SquareMatrix<double>::Zero(K, K);


    // KISS-implementation
    for (const auto& p : orbital_space.allIndices()) {
        for (const auto& q : orbital_space.allIndices()) {

            if (orbital_space.isIndexOccupied(p) && orbital_space.isIndexOccupied(q)) {  // occupied-occupied block
                const size_t i = p;
                const size_t j = q;

                double sum {0.0};
                for (const auto& a : orbital_space.virtualIndices()) {
                    sum += multipliers(j, a) * G(i, a);
                }


                if (i == j) {  // diagonal occupied part
                    Pi(i, j) += 1.0 - sum;
                } else {
                    Pi(i, j) += sum;
                }

            }  // occupied-occupied block


            else if (orbital_space.isIndexOccupied(p) && orbital_space.isIndexVirtual(q)) {  // occupied-virtual block
                const size_t i = p;
                const size_t a = q;

                double first_sum {0.0};
                for (const auto& j : orbital_space.occupiedIndices()) {
                    for (const auto& b : orbital_space.virtualIndices()) {
                        first_sum += multipliers(j, b) * G(j, b);
                    }
                }


                double second_sum {0.0};
                for (const auto& j : orbital_space.occupiedIndices()) {
                    for (const auto& b : orbital_space.virtualIndices()) {
                        if ((j != i) && (b != a)) {
                            second_sum += multipliers(j, b) * (G(i, a) * G(j, b) + G(j, a) * G(i, b));
                        }
                    }
                }

                Pi(i, a) = (1 - first_sum) * G(i, a) + second_sum;
            }  // occupied-virtual


            else if (orbital_space.isIndexVirtual(p) && orbital_space.isIndexOccupied(q)) {  // virtual-occupied block
                const size_t a = p;
                const size_t i = q;

                Pi(a, i) = multipliers(i, a);
            }


            else {  // virtual-virtual block
                const size_t a = p;
                const size_t b = q;

                double sum {0.0};
                for (const auto& i : orbital_space.occupiedIndices()) {
                    sum += multipliers(i, a) * G(i, b);
                }

                Pi(a, b) = sum;
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
TwoRDM<double> QCModel::vAP1roG::calculate2RDM(const AP1roGGeminalCoefficients& G, const BlockMatrix<double>& multipliers) {

    const size_t K = G.get_K();
    const auto orbital_space = G.orbitalSpace();

    TwoRDM<double> d {K};
    d.setZero();


    auto Delta = QCModel::vAP1roG::calculateNumber2RDM(G, multipliers);
    auto Pi = QCModel::vAP1roG::calculatePair2RDM(G, multipliers);

    // KISS-implementation
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                for (size_t s = 0; s < K; s++) {

                    if ((p == r) && (q == s)) {
                        d(p, q, r, s) += 2 * Pi(p, q);
                    }

                    if ((p == q) && (r == s) && (p != r)) {
                        d(p, q, r, s) += Delta(p, r);
                    }

                    if ((p == s) && (q == r) && (p != q)) {
                        d(p, q, r, s) -= 0.5 * Delta(p, q);
                    }
                }
            }
        }
    }  // spatial orbital loops

    return d;
}


}  // namespace GQCP
