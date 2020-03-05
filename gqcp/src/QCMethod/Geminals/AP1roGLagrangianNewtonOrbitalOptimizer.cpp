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
#include "QCMethod/Geminals/AP1roGLagrangianNewtonOrbitalOptimizer.hpp"

#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationSolver.hpp"
#include "QCMethod/Geminals/PSEnvironment.hpp"
#include "QCMethod/Geminals/vAP1roG.hpp"
#include "QCModel/Geminals/vAP1roG.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param N_P                                      the number of electron pairs
 *  @param K                                        the number of spatial orbitals
 *  @param hessian_modifier                         the modifier functor that should be used when an indefinite Hessian is encountered
 *  @param oo_convergence_threshold                 the threshold used to check for convergence on the orbital gradient
 *  @param oo_maximum_number_of_iterations          the maximum number of orbital rotation iterations that may be used to achieve convergence
 *  @param pse_convergence_threshold                the threshold used to check for convergence on the geminal coefficients
 *  @param pse_maximum_number_of_iterations         the maximum number of Newton steps that may be used to achieve convergence of the PSEs
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGLagrangianNewtonOrbitalOptimizer::AP1roGLagrangianNewtonOrbitalOptimizer(const size_t N_P, const size_t K, std::shared_ptr<BaseHessianModifier> hessian_modifier, const double oo_convergence_threshold, const size_t oo_maximum_number_of_iterations, const double pse_convergence_threshold, const size_t pse_maximum_number_of_iterations) :
    AP1roGLagrangianNewtonOrbitalOptimizer(AP1roGGeminalCoefficients(N_P, K), hessian_modifier, oo_convergence_threshold, oo_maximum_number_of_iterations, pse_convergence_threshold, pse_maximum_number_of_iterations)
{}


/**
 *  @param G                                        the initial geminal coefficients
 *  @param hessian_modifier                         the modifier functor that should be used when an indefinite Hessian is encountered
 *  @param oo_convergence_threshold                 the threshold used to check for convergence
 *  @param oo_maximum_number_of_iterations          the maximum number of iterations that may be used to achieve convergence
 *  @param pse_convergence_threshold                the threshold used to check for convergence on the geminal coefficients
 *  @param pse_maximum_number_of_iterations         the maximum number of Newton steps that may be used to achieve convergence of the PSEs
 */
AP1roGLagrangianNewtonOrbitalOptimizer::AP1roGLagrangianNewtonOrbitalOptimizer(const AP1roGGeminalCoefficients& G, std::shared_ptr<BaseHessianModifier> hessian_modifier, const double oo_convergence_threshold, const size_t oo_maximum_number_of_iterations, const double pse_convergence_threshold, const size_t pse_maximum_number_of_iterations) :
    N_P (G.get_N_P()),
    G (G),
    pse_convergence_threshold (pse_convergence_threshold),
    pse_maximum_number_of_iterations (pse_maximum_number_of_iterations),
    QCMethodNewtonOrbitalOptimizer(hessian_modifier, oo_convergence_threshold, oo_maximum_number_of_iterations)
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer
 * 
 *  In the case of this uncoupled AP1roG Lagrangian orbital optimizer, the PSEs are re-solved in every iteration using the current orbitals
 */
void AP1roGLagrangianNewtonOrbitalOptimizer::prepareDMCalculation(const SQHamiltonian<double>& sq_hamiltonian) {

    // Optimize the vAP1roG wave function model in this basis and update the results.
    auto non_linear_solver = GQCP::NonLinearEquationSolver<double>::Newton();
    auto non_linear_environment = GQCP::PSEnvironment::AP1roG(sq_hamiltonian, this->G);  // the initial guess are the current geminal coefficients
    auto linear_solver = GQCP::LinearEquationSolver<double>::HouseholderQR();

    const auto qc_structure = GQCP::QCMethod::vAP1roG(sq_hamiltonian, N_P).optimize(non_linear_solver, non_linear_environment, linear_solver);

    this->G = qc_structure.groundStateParameters().geminalCoefficients();
    this->multipliers = qc_structure.groundStateParameters().lagrangeMultipliers();
    this->E = qc_structure.groundStateEnergy();
}


/**
 *  @return the current 1-DM
 */
OneRDM<double> AP1roGLagrangianNewtonOrbitalOptimizer::calculate1RDM() const {
    return GQCP::QCModel::vAP1roG::calculate1RDM(this->G, this->multipliers);
}


/**
 *  @return the current 2-DM
 */
TwoRDM<double> AP1roGLagrangianNewtonOrbitalOptimizer::calculate2RDM() const {
    return GQCP::QCModel::vAP1roG::calculate2RDM(this->G, this->multipliers);
}


/**
 *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
 * 
 *  @param sq_hamiltonian           the current Hamiltonian
 * 
 *  @return the new full set orbital generators, including the redundant parameters
 */
OrbitalRotationGenerators AP1roGLagrangianNewtonOrbitalOptimizer::calculateNewFullOrbitalGenerators(const SQHamiltonian<double>& sq_hamiltonian) const {
    return this->calculateNewFreeOrbitalGenerators(sq_hamiltonian);  // no extra step necessary
}


}  // namespace GQCP
