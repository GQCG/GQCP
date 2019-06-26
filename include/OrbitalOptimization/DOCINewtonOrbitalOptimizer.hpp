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
#ifndef GQCP_DOCINEWTONORBITALOPTIMIZER_HPP
#define GQCP_DOCINEWTONORBITALOPTIMIZER_HPP

#include "OrbitalOptimization/QCMethodNewtonOrbitalOptimizer.hpp"

#include "HamiltonianBuilder/DOCI.hpp"
#include "Mathematical/Optimization/Eigenpair.hpp"
#include "Mathematical/Optimization/EigenproblemSolverOptions.hpp"
#include "RDM/RDMCalculator.hpp"
#include "WaveFunction/WaveFunction.hpp"

#include <memory>


namespace GQCP {


/**
 *  A class that performs gradient-and-Hessian-based orbital optimization for DOCI by sequentially
 *      - solving the DOCI eigenvalue problem
 *      - solving the Newton step to find the anti-Hermitian orbital rotation parameters
 *      - rotating the underlying spatial orbital basis
 */
class DOCINewtonOrbitalOptimizer : public QCMethodNewtonOrbitalOptimizer {
private:
    BaseSolverOptions& ci_solver_options;  // the options for the CI solver (i.e. diagonalization of the Hamiltonian)
    DOCI doci;  // the DOCI Hamiltonian builder

    RDMCalculator rdm_calculator;

    std::vector<Eigenpair> eigenpairs;  // eigenvalues and -vectors


public:
    // CONSTRUCTORS

    /**
     *  @param doci                             the DOCI HamiltonianBuilder
     *  @param ci_solver_options                the options for the CI solver (i.e. diagonalization of the Hamiltonian)
     *  @param hessian_modifier                 the modifier functor that should be used when an indefinite Hessian is encountered
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     */
    DOCINewtonOrbitalOptimizer(const DOCI& doci, BaseSolverOptions& ci_solver_options, std::shared_ptr<BaseHessianModifier> hessian_modifier, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128);


    // GETTERS

    const std::vector<Eigenpair>& get_eigenpairs() const;
    const Eigenpair& get_eigenpair(size_t index = 0) const;


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer
     * 
     *  In the case of this uncoupled DOCI orbital optimizer, the DOCI eigenvalue problem is re-solved in every iteration using the current orbitals
     */
    void prepareDMCalculation(const HamiltonianParameters<double>& ham_par) override;

    /**
     *  @return the current 1-DM
     */
    OneRDM<double> calculate1RDM() const override;

    /**
     *  @return the current 2-DM
     */
    TwoRDM<double> calculate2RDM() const override;

    /**
     *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the new full set orbital generators, including the redundant parameters
     */
    OrbitalRotationGenerators calculateNewFullOrbitalGenerators(const HamiltonianParameters<double>& ham_par) const override;


    // PUBLIC METHODS

    /**
     *  @param index        the index of the index-th excited state
     *
     *  @return the index-th excited state after doing the OO-DOCI calculation
     */
    WaveFunction makeWavefunction(size_t index = 0) const;
};


}  // namespace GQCP


#endif  // GQCP_DOCINEWTONORBITALOPTIMIZER_HPP
