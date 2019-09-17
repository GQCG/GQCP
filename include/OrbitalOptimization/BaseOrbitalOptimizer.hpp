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


#include "Basis/TransformationMatrix.hpp"
#include "Basis/SingleParticleBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  The base class for orbital optimizers. Due to the generality of the nature of the orbital optimization problem, the main algorithm (see solve()) is implemented inside this base class
 */
class BaseOrbitalOptimizer {
protected:
    bool is_converged = false;  // if the algorithm has converged
    double convergence_threshold = 1.0e-08;  // the threshold used to check for convergence
    size_t maximum_number_of_iterations = 128;  // the maximum number of iterations that may be used to achieve convergence
    size_t number_of_iterations = 0;  // the number of performed iterations


public:
    // CONSTRUCTORS

    /*
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     */
    BaseOrbitalOptimizer(const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128);


    // DESTRUCTOR
    virtual ~BaseOrbitalOptimizer() = default;


    // GETTERS
    size_t get_number_of_iterations() const { return this->number_of_iterations; }


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
     */
    virtual void prepareConvergenceChecking(const SQHamiltonian<double>& sq_hamiltonian) = 0;

    /**
     *  @param sq_hamiltonian      the current Hamiltonian
     * 
     *  @return if the algorithm is considered to be converged
     */
    virtual bool checkForConvergence(const SQHamiltonian<double>& sq_hamiltonian) const = 0;

    /**
     *  @param sq_hamiltonian      the current Hamiltonian
     * 
     *  @return a unitary matrix that will be used to rotate the current Hamiltonian into the next iteration
     */
    virtual TransformationMatrix<double> calculateNewRotationMatrix(const SQHamiltonian<double>& sq_hamiltonian) const = 0;


    // PUBLIC METHODS

    /**
     *  Optimize the Hamiltonian by subsequently
     *      - checking for convergence (see checkForConvergence())
     *      - rotating the Hamiltonian (and single-particle basis) with a newly found rotation matrix (see calculateNewRotationMatrix())
     * 
     *  @param sp_basis             the initial single-particle basis that contains the spinors to be optimized
     *  @param sq_hamiltonian       the initial (guess for the) Hamiltonian
     */
    void optimize(SingleParticleBasis<double, GTOShell>& sp_basis, SQHamiltonian<double>& sq_hamiltonian);
};


}  // namespace GQCP
