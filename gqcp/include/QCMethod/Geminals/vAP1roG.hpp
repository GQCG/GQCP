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

#pragma once


#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationEnvironment.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/QCStructure.hpp"
#include "QCModel/Geminals/vAP1roG.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  The variationally determined AP1roG quantum chemical method.
 */
class vAP1roG {
private:
    size_t N_P;  // the number of electron pairs
    size_t K;    // the number of spatial orbitals

    SQHamiltonian<double> sq_hamiltonian;  // the second-quantized Hamiltonian in an orthonormal basis


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param sq_hamiltonian           the second-quantized Hamiltonian in an orthonormal basis
     *  @param N_P                      the number of electron pairs
     */
    vAP1roG(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P) :
        N_P {N_P},
        K {sq_hamiltonian.numberOfOrbitals()},  // number of spatial orbitals
        sq_hamiltonian {sq_hamiltonian} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Optimize the vAP1roG wave function model: find the parameters satisfy the given objective.
     * 
     *  @tparam NonLinearSolver             the type of the non-linear equation solver that attempts to solve the PSEs
     *  @tparam LinearSolver                the type of the linear equation solver that attempts to determine associated Lagrange multipliers
     * 
     *  @param non_linear_solver            the solver that will try to optimize the geminal coefficients parameters
     *  @param non_linear_environment       the environment that acts as a sort of calculation space for the non-linear solver
     *  @param linear_solver                the solver that will try to find the associated Lagrange multipliers
     */
    template <typename NonLinearSolver, typename LinearSolver>
    QCStructure<GQCP::QCModel::vAP1roG, double> optimize(NonLinearSolver& non_linear_solver, NonLinearEquationEnvironment<double>& non_linear_environment, LinearSolver& linear_solver) const {

        // The vAP1roG method's responsibility is to try to:
        //      1) optimize the geminal coefficients
        //      2) find the associated Lagrange multipliers

        // 1. Optimize the geminal coefficients.
        non_linear_solver.perform(non_linear_environment);
        const auto G_optimal = AP1roGGeminalCoefficients::FromColumnMajor(non_linear_environment.variables.back(), this->N_P, this->K);


        // 2. Find the associated Lagrange multipliers.
        const auto A = QCModel::vAP1roG::calculateMultiplierResponseForceConstant(sq_hamiltonian, G_optimal);  // k_lambda
        const auto b = QCModel::vAP1roG::calculateMultiplierResponseForce(sq_hamiltonian, N_P).asVector();     // -F_lambda; column-major
        LinearEquationEnvironment<double> linear_environment {A, b};

        linear_solver.perform(linear_environment);
        const auto& lambda_optimal_vector = linear_environment.x;  // since b was column-major, so is this vector

        const size_t rows = this->N_P;            // the number of rows in the multiplier matrix
        const size_t cols = this->K - this->N_P;  // the number of columns in the multiplier matrix
        const auto orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, this->N_P}, {OccupationType::k_virtual, this->K - this->N_P}});

        const MatrixX<double> lambda_optimal_matrix = MatrixX<double>::FromColumnMajorVector(lambda_optimal_vector, rows, cols);  // the actual Lagrange multipliers, reshaped into a matrix
        const auto lambda_optimal = orbital_space.createRepresentableObjectFor<double>(OccupationType::k_occupied, OccupationType::k_virtual, lambda_optimal_matrix);


        // To make a QCStructure, we need the electronic energy, geminal coefficients and Lagrange multipliers.
        // Furthermore, the solvers only find the ground state wave function parameters, so the QCStructure only needs to contain the parameters for one state.
        const GQCP::QCModel::vAP1roG vap1rog_parameters {G_optimal, lambda_optimal};
        const auto E_electronic = vap1rog_parameters.calculateEnergy(sq_hamiltonian);

        return QCStructure<GQCP::QCModel::vAP1roG, double>({E_electronic}, {vap1rog_parameters});
    }
};


}  // namespace QCMethod
}  // namespace GQCP
