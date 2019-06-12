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
#include "OrbitalOptimization/DOCINewtonOrbitalOptimizer.hpp"

#include "CISolver/CISolver.hpp"
#include "math/optimization/step.hpp"
#include "RDM/RDMCalculator.hpp"
#include "utilities/linalg.hpp"



namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param doci                     the DOCI HamiltonianBuilder
 *  @param ci_solver_options        the options for the CI solver (i.e. diagonalization of the Hamiltonian)
 *  @param oo_options               the options for orbital optimization
 */
DOCINewtonOrbitalOptimizer::DOCINewtonOrbitalOptimizer(const DOCI& doci, BaseSolverOptions& ci_solver_options, const OrbitalOptimizationOptions& oo_options) :
    doci (doci),
    ci_solver_options (ci_solver_options),
    rdm_calculator (RDMCalculator(*this->doci.get_fock_space())),
    NewtonOrbitalOptimizer(oo_options)
{}


/*
 *  GETTERS
 */

const std::vector<Eigenpair>& DOCINewtonOrbitalOptimizer::get_eigenpairs() const {
    if (this->is_converged) {
        return this->eigenpairs;
    } else {
        throw std::logic_error("DOCINewtonOrbitalOptimizer::get_eigenpairs(): You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
    }
}

const Eigenpair& DOCINewtonOrbitalOptimizer::get_eigenpair(size_t index) const {
    if (this->is_converged) {
        return this->eigenpairs[index];
    } else {
        throw std::logic_error("DOCINewtonOrbitalOptimizer::get_eigenpair(size_t): You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
    }
}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */


/**
 *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer
 * 
 *  In the case of this uncoupled DOCI orbital optimizer, the DOCI eigenvalue problem is re-solved in every iteration using the current orbitals
 */
void DOCINewtonOrbitalOptimizer::prepareNewtonSpecificConvergenceChecking(const HamiltonianParameters<double>& ham_par) {

    // Solve the DOCI eigenvalue problem to obtain DMs from which we can calculate the gradient and the Hessian
    CISolver doci_solver (this->doci, ham_par);
    doci_solver.solve(this->ci_solver_options);
    this->eigenpairs = doci_solver.get_eigenpairs();

    // If we're using a Davidson solver, we should update the initial guesses in the CI solver options to be the current eigenvectors
    if (this->ci_solver_options.get_solver_type() == SolverType::DAVIDSON) {
        auto davidson_solver_options = dynamic_cast<DavidsonSolverOptions&>(this->ci_solver_options);

        for (size_t i = 0; i < davidson_solver_options.number_of_requested_eigenpairs; i++) {
            davidson_solver_options.X_0.col(i) = doci_solver.makeWavefunction(i).get_coefficients();
        }

        this->ci_solver_options = davidson_solver_options;
    }


    // Calculate the 1- and 2-DMs
    this->rdm_calculator.set_coefficients(doci_solver.get_eigenpair().get_eigenvector());
    this->D = this->rdm_calculator.calculate1RDMs().one_rdm;
    this->d = this->rdm_calculator.calculate2RDMs().two_rdm;
}


/**
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return the current orbital gradient as a matrix
 */
SquareMatrix<double> DOCINewtonOrbitalOptimizer::calculateGradientMatrix(const HamiltonianParameters<double>& ham_par) const {

    // Calculate the gradient from the Fockian matrix
    const auto F = ham_par.calculateGeneralizedFockMatrix(this->D, this->d);
    return 2 * (F - F.transpose());
}


/**
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return the current orbital Hessian as a tensor
 */
SquareRankFourTensor<double> DOCINewtonOrbitalOptimizer::calculateHessianTensor(const HamiltonianParameters<double>& ham_par) const {

    const auto K = ham_par.get_K();

    // The 1- and 2-DM have already been calculated in this->calculateGradientMatrix
    const auto G = ham_par.calculateSuperGeneralizedFockMatrix(this->D, this->d);
    SquareRankFourTensor<double> hessian_tensor (K);
    hessian_tensor.setZero();

    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                for (size_t s = 0; s < K; s++) {
                    hessian_tensor(p,q,r,s) = G(p,q,r,s) - G(p,q,s,r) + G(q,p,s,r) - G(q,p,r,s) + G(r,s,p,q) - G(r,s,q,p) + G(s,r,q,p) - G(s,r,p,q);
                }
            }
        }
    }

    return hessian_tensor;
}


/**
 *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
 * 
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return the new full set orbital generators, including the redundant parameters
 */
OrbitalRotationGenerators DOCINewtonOrbitalOptimizer::calculateNewFullOrbitalGenerators(const HamiltonianParameters<double>& ham_par) const {
    return this->calculateNewFreeOrbitalGenerators(ham_par);  // no extra step necessary
}




/*
 *  PUBLIC METHODS
 */

/**
 *  @param index        the index of the index-th excited state
 *
 *  @return the index-th excited state after doing the OO-DOCI calculation
 */
WaveFunction DOCINewtonOrbitalOptimizer::makeWavefunction(size_t index) const {

    if (index > this->eigenpairs.size()) {
        throw std::logic_error("DOCINewtonOrbitalOptimizer::makeWavefunction(size_t): Not enough requested eigenpairs for the given index.");
    }

    return WaveFunction(*this->doci.get_fock_space(), this->eigenpairs[index].get_eigenvector());
}


}  // namespace GQCP
