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


#include "HamiltonianBuilder/DOCI.hpp"
#include "math/optimization/Eigenpair.hpp"
#include "math/optimization/EigenproblemSolverOptions.hpp"
#include "OrbitalOptimization/NewtonOrbitalOptimizer.hpp"
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
class DOCINewtonOrbitalOptimizer : public NewtonOrbitalOptimizer {
private:
    BaseSolverOptions& ci_solver_options;  // the options for the CI solver (i.e. diagonalization of the Hamiltonian)
    DOCI doci;  // the DOCI Hamiltonian builder

    RDMCalculator rdm_calculator;
    OneRDM<double> D;  // spin-summed 1-RDM
    TwoRDM<double> d;  // spin-summed 2-RDM

    std::vector<Eigenpair> eigenpairs;  // eigenvalues and -vectors


public:
    // CONSTRUCTORS

    /**
     *  @param doci                     the DOCI HamiltonianBuilder
     *  @param ci_solver_options        the options for the CI solver (i.e. diagonalization of the Hamiltonian)
     *  @param oo_options               the options for orbital optimization
     */
    DOCINewtonOrbitalOptimizer(const DOCI& doci, BaseSolverOptions& ci_solver_options, const OrbitalOptimizationOptions& oo_options);


    // GETTERS

    const std::vector<Eigenpair>& get_eigenpairs() const;
    const Eigenpair& get_eigenpair(size_t index = 0) const;


    // OVERRIDDEN PUBLIC METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer
     * 
     *  In the case of this uncoupled DOCI orbital optimizer, the DOCI eigenvalue problem is re-solved in every iteration using the current orbitals
     */
    void prepareNewtonSpecificConvergenceChecking(const HamiltonianParameters<double>& ham_par) override;

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to calculate the new rotation matrix in this Newton-based orbital optimizer
     */
    void prepareNewtonSpecificRotationMatrixCalculation(const HamiltonianParameters<double>& ham_par) override {}

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital gradient as a matrix
     */
    SquareMatrix<double> calculateGradientMatrix(const HamiltonianParameters<double>& ham_par) const override;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital Hessian as a tensor
     */
    SquareRankFourTensor<double> calculateHessianTensor(const HamiltonianParameters<double>& ham_par) const override;

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
