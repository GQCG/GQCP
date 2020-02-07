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
#include "QCMethod/Geminals/AP1roGPSESolver.hpp"

#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationSolver.hpp"
#include "QCMethod/Geminals/AP1roG.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  @param pses         the AP1roG PSEs
 */
AP1roGPSESolver::AP1roGPSESolver(const AP1roGPSEs& pses, const double convergence_threshold, const size_t maximum_number_of_iterations) : 
    convergence_threshold (convergence_threshold),
    maximum_number_of_iterations (maximum_number_of_iterations),
    pses (pses)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the projected Schrödinger equations for AP1roG
 * 
 *  @param G            the initial geminal coefficients, that are updated in every iteration to the converged geminal coefficients
 */
void AP1roGPSESolver::solve(AP1roGGeminalCoefficients& G) const {

    const auto f = this->pses.callableCoordinateFunctions();
    const auto J = this->pses.callableJacobian();

    // Solve the AP1roG equations using a Newton-based algorithm
    VectorX<double> x = G.asVector();  // the initial guess: a column-major vector
    NonLinearEquationEnvironment<double> non_linear_environment (x, f, J);
    auto non_linear_solver = NonLinearEquationSolver<double>::Newton(this->convergence_threshold, this->maximum_number_of_iterations);
    non_linear_solver.perform(non_linear_environment);
    x = non_linear_environment.variables.back();

    G = AP1roGGeminalCoefficients::FromColumnMajor(x, this->pses.numberOfElectronPairs(), this->pses.numberOfSpatialOrbitals());
}



/**
 *  Solve the projected Schrödinger equations for AP1roG, using a zero initial guess
 */
AP1roGGeminalCoefficients AP1roGPSESolver::solve() const {

    const auto N_P = this->pses.numberOfElectronPairs();
    const auto K = this->pses.numberOfSpatialOrbitals();

    // Set up the zero initial guess and pass it to the other solver method that accepts an initial guess
    AP1roGGeminalCoefficients G (N_P, K);
    this->solve(G);
    return G;
}

}  // namespace GQCP
