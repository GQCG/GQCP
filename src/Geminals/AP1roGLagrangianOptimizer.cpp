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

#include "Geminals/AP1roG.hpp"
#include "Geminals/AP1roGPSEs.hpp"
#include "Geminals/AP1roGPSESolver.hpp"


namespace GQCP {


/*
 *  PUBLIC METHODS
 */
void AP1roGLagrangianOptimizer::solve() {

    const auto& g = this->ham_par.get_g();
    auto K = this->ham_par.get_K();


    // Solve the PSEs and set part of the solutions
    AP1roGPSEs pses (this->ham_par, this->N_P);
    AP1roGPSESolver pse_solver (pses);
    pse_solver.solve(this->geminal_coefficients);  // changes the geminal coefficients into the solution

    this->electronic_energy = calculateAP1roGEnergy(this->geminal_coefficients, this->ham_par);


    // Initialize and solve the linear system k_lambda lambda=b (lambda are the Lagrange multipliers)
    // b_{ia} = dE/dG_i^a
    const MatrixX<double> k_lambda = pses.evaluateJacobian(this->geminal_coefficients).asMatrix().transpose();

    BlockMatrix<double> b (0, this->N_P, this->N_P, K);
    for (size_t i = 0; i < this->N_P; i++) {
        for (size_t a = this->N_P; a < this->K; a++) {
            b(i,a) = g(i,a,i,a);
        }
    }

    Eigen::HouseholderQR<Eigen::MatrixXd> linear_solver (k_lambda);
    VectorX<double> lambda = linear_solver.solve(-b.asVector());


    if (std::abs((k_lambda * lambda + b.asVector()).norm()) > 1.0e-12) {
        throw std::runtime_error("void AP1roGLagrangianOptimizer::solve(): The Householder QR decomposition failed.");
    }

    this->multipliers = AP1roGVariables(lambda, this->N_P, K);
}


}  // namespace GQCP
