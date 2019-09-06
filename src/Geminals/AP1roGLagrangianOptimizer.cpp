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
#include "Geminals/AP1roGLagrangianOptimizer.hpp"

#include "Geminals/AP1roGPSESolver.hpp"


namespace GQCP {


/*
 *  PUBLIC METHODS
 */
void AP1roGLagrangianOptimizer::solve() {

    const auto K = this->ham_par.get_K();
    const auto& g = this->ham_par.twoElectron().parameters();


    // Solve the PSEs and set part of the solutions
    AP1roGPSESolver pse_solver (this->N_P, this->ham_par, this->geminal_coefficients, this->convergence_threshold, this->maximum_number_of_iterations);
    pse_solver.solve();

    this->geminal_coefficients = pse_solver.get_geminal_coefficients();
    this->electronic_energy = pse_solver.get_electronic_energy();


    // Initialize and solve the linear system Jx=b (x are the Lagrange multipliers)
    size_t dim = this->geminal_coefficients.numberOfGeminalCoefficients(this->N_P, K);

    auto J = pse_solver.calculateJacobian(this->geminal_coefficients);

    VectorX<double> b = VectorX<double>::Zero(dim);  // dE/dG_i^a
    for (size_t i = 0; i < this->N_P; i++) {
        for (size_t a = this->N_P; a < this->K; a++) {
            size_t row_vector_index = this->geminal_coefficients.vectorIndex(i, a);

            b(row_vector_index) = -g(i,a,i,a);
        }
    }

    Eigen::HouseholderQR<Eigen::MatrixXd> linear_solver (J);
    VectorX<double> x = linear_solver.solve(b);

    if (std::abs((J * x - b).norm()) > 1.0e-12) {
        throw std::runtime_error("void AP1roGLagrangianOptimizer::solve(): The Householder QR decomposition failed.");
    }

    this->multipliers = AP1roGVariables(x, this->N_P, K);
}


}  // namespace GQCP
