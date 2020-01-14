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
#include "QCMethod/Geminals/AP1roGLagrangianOptimizer.hpp"

#include "QCMethod/Geminals/AP1roGPSEs.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param G                    the converged geminal coefficients that are a solution to the AP1roG PSEs
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal orbital basis
 */
AP1roGLagrangianOptimizer::AP1roGLagrangianOptimizer(const AP1roGGeminalCoefficients& G, const SQHamiltonian<double>& sq_hamiltonian) :
    G (G),
    sq_hamiltonian (sq_hamiltonian)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the Lagrange multipliers for the AP1roG PSE Lagrangian
 */
BlockMatrix<double> AP1roGLagrangianOptimizer::solve() {

    // Initialize some variables that are needed in the solution algorithm
    const auto K = this->G.numberOfSpatialOrbitals();
    const auto N_P = this->G.numberOfElectronPairs();
    const auto& g = this->sq_hamiltonian.twoElectron().parameters();


    // Initialize and solve the linear system (k_lambda lambda=-b) in order to determine the Lagrange multipliers lambda

    // Calculate the matrix k_lambda
    const AP1roGPSEs pses (this->sq_hamiltonian, N_P);
    const MatrixX<double> k_lambda = pses.calculateJacobian(this->G).asMatrix().transpose();


    // Calculate the right-hand side vector b_i^a = dE/DG_i^a
    BlockMatrix<double> b (0, N_P, N_P, K);
    for (size_t i = 0; i < N_P; i++) {
        for (size_t a = N_P; a < K; a++) {
            b(i,a) = g(i,a,i,a);
        }
    }


    // Solve the linear system (k_lambda lambda=-b)
    Eigen::HouseholderQR<Eigen::MatrixXd> linear_solver (k_lambda);
    VectorX<double> lambda = linear_solver.solve(-b.asVector());  // b is column major, so lambda is also column major


    if (std::abs((k_lambda * lambda + b.asVector()).norm()) > 1.0e-12) {
        throw std::runtime_error("void AP1roGLagrangianOptimizer::solve(): The Householder QR decomposition failed.");
    }

    const size_t rows = N_P;  // the number of rows in the multiplier matrix
    const size_t cols = K - N_P;  // the number of columns in the multiplier matrix
    const MatrixX<double> M = MatrixX<double>::FromColumnMajorVector(lambda, rows, cols);  // the actual Lagrange multipliers, reshaped into a matrix
    return BlockMatrix<double>(0, N_P, N_P, K, M);  // an encapsulating object that implements operator() intuitively
}


}  // namespace GQCP
