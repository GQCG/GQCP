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
#include "Mathematical/Optimization/NewtonSystemOfEquationsSolver.hpp"


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
 *  @param G            the initial geminal coefficients
 */
AP1roGGeminalCoefficients AP1roGPSESolver::solve(const AP1roGGeminalCoefficients& G) const {

    // Solve the AP1roG equations using a Newton-based algorithm
    const auto f = this->pses.callableCoordinateFunctions();
    const auto J = this->pses.callableJacobian();


    VectorX<double> x0 = G.asVector();
    NewtonSystemOfEquationsSolver syseq_solver (x0, f, J);
    syseq_solver.solve();

    const auto N_P = this->pses.numberOfElectronPairs();
    const auto K = this->pses.numberOfSpatialOrbitals();

    return AP1roGGeminalCoefficients::FromColumnMajor(syseq_solver.get_solution(), N_P, K);
}


/**
 *  Solve the projected Schrödinger equations for AP1roG, using a zero initial guess
 */
AP1roGGeminalCoefficients AP1roGPSESolver::solve() const {

    const auto N_P = this->pses.numberOfElectronPairs();
    const auto K = this->pses.numberOfSpatialOrbitals();
    AP1roGGeminalCoefficients G_initial (N_P, K);

    return this->solve(G_initial);
}


}  // namespace GQCP
