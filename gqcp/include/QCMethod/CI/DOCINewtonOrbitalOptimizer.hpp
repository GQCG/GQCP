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
#pragma once


#include "Mathematical/Optimization/Eigenproblem/Eigenpair.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "Processing/RDM/DOCIRDMBuilder.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/OrbitalOptimization/QCMethodNewtonOrbitalOptimizer.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <memory>


namespace GQCP {


/**
 *  A class that performs gradient-and-Hessian-based orbital optimization for DOCI by sequentially
 *      - solving the DOCI eigenvalue problem
 *      - solving the Newton step to find the anti-Hermitian orbital rotation parameters
 *      - rotating the underlying spatial orbital basis
 * 
 *  @tparam _EigenproblemSolver          the type of the eigenproblem solver that is used
 */
template <typename _EigenproblemSolver>
class DOCINewtonOrbitalOptimizer: 
    public QCMethodNewtonOrbitalOptimizer {

public:
    using EigenproblemSolver = _EigenproblemSolver;


private:
    SeniorityZeroONVBasis onv_basis;  // the Fock subspace used for DOCI calculations

    EigenproblemEnvironment eigenproblem_environment;
    EigenproblemSolver eigenproblem_solver;

    DOCIRDMBuilder rdm_calculator;

    LinearExpansion<SeniorityZeroONVBasis> ground_state_expansion;

    size_t number_of_requested_eigenpairs;
    std::vector<Eigenpair> eigenpairs;  // eigenvalues and -vectors


public:
    // CONSTRUCTORS

    /**
     *  @param onv_basis                        the Fock subspace used for DOCI calculations
     *  @param eigenproblem_environment         the environments that acts as the calculation context for the eigenproblem solver
     *  @param eigenproblem_solver              the algorithm that tries to solve the DOCI eigenvalue problem
     *  @param hessian_modifier                 the modifier functor that should be used when an indefinite Hessian is encountered
     *  @param number_of_requested_eigenpairs   the number of eigenpairs that should be looked for
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     */
    DOCINewtonOrbitalOptimizer(const SeniorityZeroONVBasis& onv_basis, const EigenproblemSolver& eigenproblem_solver, const EigenproblemEnvironment& eigenproblem_environment, std::shared_ptr<BaseHessianModifier> hessian_modifier, const size_t number_of_requested_eigenpairs = 1, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) :
        onv_basis (onv_basis),
        eigenproblem_solver (eigenproblem_solver),
        eigenproblem_environment (eigenproblem_environment),
        number_of_requested_eigenpairs (number_of_requested_eigenpairs),
        rdm_calculator (DOCIRDMBuilder(onv_basis)),
        QCMethodNewtonOrbitalOptimizer(hessian_modifier, convergence_threshold, maximum_number_of_iterations)
    {}


    // GETTERS

    const std::vector<Eigenpair>& get_eigenpairs() const {
        if (this->is_converged) {
            return this->eigenpairs;
        } else {
            throw std::logic_error("DOCINewtonOrbitalOptimizer::get_eigenpairs(): You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
        }
    }


    const Eigenpair& get_eigenpair(size_t index = 0) const {
        if (this->is_converged) {
            return this->eigenpairs[index];
        } else {
            throw std::logic_error("DOCINewtonOrbitalOptimizer::get_eigenpair(size_t): You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
        }
    }


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer
     * 
     *  In the case of this uncoupled DOCI orbital optimizer, the DOCI eigenvalue problem is re-solved in every iteration using the current orbitals
     */
    void prepareDMCalculation(const SQHamiltonian<double>& sq_hamiltonian) override {

        // (Re)create the eigenproblem environment in the current orbital basis.
        if (this->eigenproblem_environment.A.cols() != 0) {  // if the optimization environment is 'dense'
            this->eigenproblem_environment = CIEnvironment::Dense(sq_hamiltonian, this->onv_basis);
        } else {

            // Recreate the iterative eigenproblem environment with the previous guesses.
            if (this->number_of_iterations != 0) {  // not needed when we haven't done an OO iteration
                const auto dim = this->onv_basis.dimension();
                MatrixX<double> V = MatrixX<double>::Zero(dim, this->number_of_requested_eigenpairs);

                const auto eigenpairs = this->eigenproblem_environment.eigenpairs(this->number_of_requested_eigenpairs);
                for (size_t i = 0; i < this->number_of_requested_eigenpairs; i++) {
                    V.col(i) = eigenpairs[i].get_eigenvector();
                }

                this->eigenproblem_environment = CIEnvironment::Iterative(sq_hamiltonian, this->onv_basis, V);
            }
        }

        // Set the ground state expansion and the possibly requested excited states.
        this->ground_state_expansion = QCMethod::CI<SeniorityZeroONVBasis>(this->onv_basis, this->number_of_requested_eigenpairs).optimize(this->eigenproblem_solver, this->eigenproblem_environment).groundStateParameters();
        this->eigenpairs = eigenproblem_environment.eigenpairs(this->number_of_requested_eigenpairs);
    }

    /**
     *  @return the current 1-DM
     */
    OneRDM<double> calculate1RDM() const override {
        return this->rdm_calculator.calculate1RDMs(this->ground_state_expansion.coefficients()).one_rdm;
    }

    /**
     *  @return the current 2-DM
     */
    TwoRDM<double> calculate2RDM() const override {
        return this->rdm_calculator.calculate2RDMs(this->ground_state_expansion.coefficients()).two_rdm;
    }

    /**
     *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
     * 
     *  @param sq_hamiltonian      the current Hamiltonian
     * 
     *  @return the new full set orbital generators, including the redundant parameters
     */
    OrbitalRotationGenerators calculateNewFullOrbitalGenerators(const SQHamiltonian<double>& sq_hamiltonian) const override {
        return this->calculateNewFreeOrbitalGenerators(sq_hamiltonian);  // no extra step necessary
    }


    // PUBLIC METHODS

    /**
     *  @param index        the index of the index-th excited state
     *
     *  @return the index-th excited state after doing the OO-DOCI calculation
     */
    LinearExpansion<SeniorityZeroONVBasis> makeLinearExpansion(size_t index = 0) const {
        if (index > this->eigenpairs.size()) {
            throw std::logic_error("DOCINewtonOrbitalOptimizer::makeLinearExpansion(size_t): Not enough requested eigenpairs for the given index.");
        }

        return LinearExpansion<SeniorityZeroONVBasis>{this->onv_basis, this->eigenpairs[index].get_eigenvector()};
    }
};


}  // namespace GQCP
