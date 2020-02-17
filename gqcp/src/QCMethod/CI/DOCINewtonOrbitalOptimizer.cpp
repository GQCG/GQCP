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
#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"

#include "Mathematical/Optimization/NonLinearEquation/step.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "Utilities/linalg.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param doci                             the DOCI HamiltonianBuilder
 *  @param ci_solver_options                the options for the CI solver (i.e. diagonalization of the Hamiltonian)
 *  @param hessian_modifier                 the modifier functor that should be used when an indefinite Hessian is encountered
 *  @param convergence_threshold            the threshold used to check for convergence
 *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
 */
DOCINewtonOrbitalOptimizer::DOCINewtonOrbitalOptimizer(const DOCI& doci, BaseSolverOptions& ci_solver_options, std::shared_ptr<BaseHessianModifier> hessian_modifier, const double convergence_threshold, const size_t maximum_number_of_iterations) :
    doci (doci),
    ci_solver_options (ci_solver_options),
    rdm_calculator (RDMCalculator(*this->doci.get_fock_space())),
    QCMethodNewtonOrbitalOptimizer(hessian_modifier, convergence_threshold, maximum_number_of_iterations)
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
void DOCINewtonOrbitalOptimizer::prepareDMCalculation(const SQHamiltonian<double>& sq_hamiltonian) {

    // Solve the DOCI eigenvalue problem to obtain DMs from which we can calculate the gradient and the Hessian
    CISolver doci_solver (this->doci, sq_hamiltonian);
    doci_solver.solve(this->ci_solver_options);
    this->eigenpairs = doci_solver.get_eigenpairs();

    // If we're using a Davidson solver, we should update the initial guesses in the CI solver options to be the current eigenvectors
    if (this->ci_solver_options.get_solver_type() == SolverType::DAVIDSON) {
        auto davidson_solver_options = dynamic_cast<DavidsonSolverOptions&>(this->ci_solver_options);

        for (size_t i = 0; i < davidson_solver_options.number_of_requested_eigenpairs; i++) {
            davidson_solver_options.X_0.col(i) = doci_solver.makeLinearExpansion(i).get_coefficients();
        }

        this->ci_solver_options = davidson_solver_options;
    }

    // Prepare the calculation of the 1- and 2-DMs by putting in the most current coefficients
    this->rdm_calculator.set_coefficients(doci_solver.get_eigenpair().get_eigenvector());
}


/**
 *  @return the current 1-DM
 */
OneRDM<double> DOCINewtonOrbitalOptimizer::calculate1RDM() const {
    return this->rdm_calculator.calculate1RDMs().one_rdm;
}


/**
 *  @return the current 2-DM
 */
TwoRDM<double> DOCINewtonOrbitalOptimizer::calculate2RDM() const {
    return this->rdm_calculator.calculate2RDMs().two_rdm;
}


/**
 *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
 * 
 *  @param sq_hamiltonian           the current Hamiltonian
 * 
 *  @return the new full set orbital generators, including the redundant parameters
 */
OrbitalRotationGenerators DOCINewtonOrbitalOptimizer::calculateNewFullOrbitalGenerators(const SQHamiltonian<double>& sq_hamiltonian) const {
    return this->calculateNewFreeOrbitalGenerators(sq_hamiltonian);  // no extra step necessary
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param index        the index of the index-th excited state
 *
 *  @return the index-th excited state after doing the OO-DOCI calculation
 */
LinearExpansion DOCINewtonOrbitalOptimizer::makeLinearExpansion(size_t index) const {

    if (index > this->eigenpairs.size()) {
        throw std::logic_error("DOCINewtonOrbitalOptimizer::makeLinearExpansion(size_t): Not enough requested eigenpairs for the given index.");
    }

    return LinearExpansion(*this->doci.get_fock_space(), this->eigenpairs[index].get_eigenvector());
}


}  // namespace GQCP
