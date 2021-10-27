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


#include "Basis/Transformations/ROrbitalRotationGenerators.hpp"
#include "Mathematical/Optimization/Eigenproblem/Eigenpair.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/OrbitalOptimization/QCMethodNewtonOrbitalOptimizer.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <memory>


namespace GQCP {


/**
 * MARK: DOCINewtonOrbitalOPtimizer implementation
 */

/**
 *  A class that performs gradient-and-Hessian-based orbital optimization for DOCI by sequentially:
 *      - solving the DOCI eigenvalue problem,
 *      - solving the Newton step to find the anti-Hermitian orbital rotation parameters,
 *      - rotating the underlying spatial orbital basis.
 *
 *  @tparam _EigenproblemSolver          The type of the eigenproblem solver that is used.
 */
template <typename _EigenproblemSolver>
class DOCINewtonOrbitalOptimizer:
    public QCMethodNewtonOrbitalOptimizer {
public:
    // The type of Eigenproblem solver.
    using EigenproblemSolver = _EigenproblemSolver;

    // The type of Eigenproblem environment.
    using EigenproblemEnvironment = EigenproblemEnvironment<double>;

    // The type of eigenpairs used in this orbital optimizer.
    using Eigenpair = Eigenpair<double, double>;

    // The type of ONV basis used for DOCI calculations.
    using ONVBasis = SeniorityZeroONVBasis;

    // The type of linear expansion wave function model.
    using LinearExpansion = LinearExpansion<double, ONVBasis>;


private:
    // The Fock subspace used for DOCI calculations.
    ONVBasis onv_basis;

    // The associated eigenproblem solver and environment.
    EigenproblemEnvironment eigenproblem_environment;
    EigenproblemSolver eigenproblem_solver;

    // The linear expansion representing the DOCI wave function model.
    LinearExpansion ground_state_expansion;

    // The number of requested eigenpairs.
    size_t number_of_requested_eigenpairs;

    // The eigenpairs containing the eigenvalues and eigenvectors of the eigenvalue problem.
    std::vector<Eigenpair> m_eigenpairs;


public:
    /**
     * MARK: Constructors
     */

    /**
     *  @param onv_basis                        The Fock subspace used for DOCI calculations.
     *  @param eigenproblem_environment         The environments that acts as the calculation context for the eigenproblem solver.
     *  @param eigenproblem_solver              The algorithm that tries to solve the DOCI eigenvalue problem.
     *  @param hessian_modifier                 The modifier functor that should be used when an indefinite Hessian is encountered.
     *  @param number_of_requested_eigenpairs   The number of m_eigenpairs that should be looked for.
     *  @param convergence_threshold            The threshold used to check for convergence.
     *  @param maximum_number_of_iterations     The maximum number of iterations that may be used to achieve convergence.
     */
    DOCINewtonOrbitalOptimizer(const SeniorityZeroONVBasis& onv_basis, const EigenproblemSolver& eigenproblem_solver, const EigenproblemEnvironment& eigenproblem_environment, std::shared_ptr<BaseHessianModifier> hessian_modifier, const size_t number_of_requested_eigenpairs = 1, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) :
        onv_basis {onv_basis},
        eigenproblem_solver {eigenproblem_solver},
        eigenproblem_environment {eigenproblem_environment},
        number_of_requested_eigenpairs {number_of_requested_eigenpairs},
        QCMethodNewtonOrbitalOptimizer(hessian_modifier, convergence_threshold, maximum_number_of_iterations) {}


    /**
     * MARK: Overridden methods
     */

    /**
     *  @return The current 1-DM.
     */
    Orbital1DM<double> calculate1DM() const override {
        return this->ground_state_expansion.calculate1DM();
    }


    /**
     *  @return The current 2-DM.
     */
    Orbital2DM<double> calculate2DM() const override {
        return this->ground_state_expansion.calculate2DM();
    }


    /**
     *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential.
     *
     *  @param sq_hamiltonian      The current second quantized Hamiltonian.
     *
     *  @return The new full set orbital generators, including the redundant parameters.
     */
    ROrbitalRotationGenerators<double> calculateNewFullOrbitalGenerators(const RSQHamiltonian<double>& sq_hamiltonian) const override {
        return this->calculateNewFreeOrbitalGenerators(sq_hamiltonian);  // No extra step necessary.
    }


    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer.
     *
     *  In the case of this uncoupled DOCI orbital optimizer, the DOCI eigenvalue problem is re-solved in every iteration using the current orbitals.
     */
    void prepareDMCalculation(const RSQHamiltonian<double>& sq_hamiltonian) override {

        // (Re)create the eigenproblem environment in the current orbital basis.
        if (this->eigenproblem_environment.A.cols() != 0) {  // If the optimization environment is 'dense'.
            this->eigenproblem_environment = CIEnvironment::Dense(sq_hamiltonian, this->onv_basis);
        } else {

            // Recreate the iterative eigenproblem environment with the previous guesses.
            if (this->number_of_iterations != 0) {  // Not needed when we haven't done an OO iteration.
                const auto dim = this->onv_basis.dimension();
                MatrixX<double> V = MatrixX<double>::Zero(dim, this->number_of_requested_eigenpairs);

                const auto m_eigenpairs = this->eigenproblem_environment.eigenpairs(this->number_of_requested_eigenpairs);
                for (size_t i = 0; i < this->number_of_requested_eigenpairs; i++) {
                    V.col(i) = m_eigenpairs[i].eigenvector();
                }

                this->eigenproblem_environment = CIEnvironment::Iterative(sq_hamiltonian, this->onv_basis, V);
            }
        }

        // Set the ground state expansion and the possibly requested excited states.
        this->ground_state_expansion = QCMethod::CI<double, ONVBasis>(this->onv_basis, this->number_of_requested_eigenpairs).optimize(this->eigenproblem_solver, this->eigenproblem_environment).groundStateParameters();
        this->m_eigenpairs = eigenproblem_environment.eigenpairs(this->number_of_requested_eigenpairs);
    }


    /**
     * MARK: Eigenpair access
     */

    /**
     *  @param index                The index of a state.
     *
     *  @return The eigenpair that is associated to the given index.
     */
    const Eigenpair& eigenpair(const size_t index = 0) const {

        if (this->is_converged) {
            return this->m_eigenpairs[index];
        } else {
            throw std::logic_error("DOCINewtonOrbitalOptimizer::eigenpair(const size_t): You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
        }
    }

    /**
     *  @return All eigenpairs found by this orbital optimizer.
     */
    const std::vector<Eigenpair>& eigenpairs() const {

        if (this->is_converged) {
            return this->m_eigenpairs;
        } else {
            throw std::logic_error("DOCINewtonOrbitalOptimizer::eigenpairs(): You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
        }
    }

    /**
     * MARK: Wave function model
     */

    /**
     *  @param index        The index of the index-th excited state.
     *
     *  @return The index-th excited state after doing the OO-DOCI calculation.
     */
    LinearExpansion makeLinearExpansion(size_t index = 0) const {
        if (index > this->m_eigenpairs.size()) {
            throw std::logic_error("DOCINewtonOrbitalOptimizer::makeLinearExpansion(size_t): Not enough requested m_eigenpairs for the given index.");
        }

        return LinearExpansion(this->onv_basis, this->m_eigenpairs[index].eigenvector());
    }
};


}  // namespace GQCP
