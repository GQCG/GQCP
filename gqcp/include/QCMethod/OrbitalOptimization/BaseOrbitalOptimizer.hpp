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


#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  The base class for orbital optimizers. Due to the generality of the nature of the orbital optimization problem, the main algorithm (see solve()) is implemented inside this base class
 */
class BaseOrbitalOptimizer {
protected:
    bool is_converged = false;                  // if the algorithm has converged
    double convergence_threshold = 1.0e-08;     // the threshold used to check for convergence
    size_t maximum_number_of_iterations = 128;  // the maximum number of iterations that may be used to achieve convergence
    size_t number_of_iterations = 0;            // the number of performed iterations


public:
    // CONSTRUCTORS

    /*
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     */
    BaseOrbitalOptimizer(const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128);


    // DESTRUCTOR
    virtual ~BaseOrbitalOptimizer() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  @param sq_hamiltonian      the current Hamiltonian
     * 
     *  @return The unitary transformation that will be used to rotate the current Hamiltonian into the next iteration.
     */
    virtual RTransformation<double> calculateNewRotationMatrix(const RSQHamiltonian<double>& sq_hamiltonian) const = 0;

    /**
     *  @param sq_hamiltonian      the current Hamiltonian
     * 
     *  @return if the algorithm is considered to be converged
     */
    virtual bool checkForConvergence(const RSQHamiltonian<double>& sq_hamiltonian) const = 0;

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
     */
    virtual void prepareConvergenceChecking(const RSQHamiltonian<double>& sq_hamiltonian) = 0;


    // PUBLIC METHODS

    /**
     *  @return the number of iterations that this optimizer has performed
     */
    size_t numberOfIterations() const { return this->number_of_iterations; }

    /**
     *  Optimize the Hamiltonian by subsequently
     *      - checking for convergence (see checkForConvergence())
     *      - rotating the Hamiltonian (and spinor basis) with a newly found rotation matrix (see calculateNewRotationMatrix())
     * 
     *  @param spinor_basis         the initial spinor basis that contains the spinors to be optimized
     *  @param sq_hamiltonian       the initial (guess for the) Hamiltonian
     */
    void optimize(RSpinOrbitalBasis<double, GTOShell>& spinor_basis, RSQHamiltonian<double>& sq_hamiltonian);
};


}  // namespace GQCP
