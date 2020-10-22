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


#include "Basis/Transformations/JacobiRotationParameters.hpp"
#include "QCMethod/OrbitalOptimization/BaseOrbitalOptimizer.hpp"


namespace GQCP {


/**
 *  An intermediate abstract class that should be derived from to implement a Jacobi rotation based orbital optimization: the change in scalar function due to a Jacobi rotation should be implemented
 */
class JacobiOrbitalOptimizer: public BaseOrbitalOptimizer {
protected:
    size_t dim;  // the dimension of the orbital space that should be scanned. The valid orbital indices then are 0 ... dim (not included)

    using pair_type = std::pair<JacobiRotationParameters, double>;
    pair_type optimal_jacobi_with_scalar;  // holds the optimal Jacobi parameters and the corresponding value for the scalar function trying to optimize


public:
    // CONSTRUCTORS

    /**
     *  @param dim                             the dimension of the orbital space that should be scanned. The valid orbital indices then are 0 ... dim (not included)
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     */
    JacobiOrbitalOptimizer(const size_t dim, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128);


    // DESTRUCTOR

    /**
     *  The default destructor.
     */
    virtual ~JacobiOrbitalOptimizer() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  Calculate the trigoniometric polynomial coefficients for the given Jacobi rotation indices
     *
     *  @param p            the index of spatial orbital 1
     *  @param q            the index of spatial orbital 2
     */
    virtual void calculateJacobiCoefficients(const RSQHamiltonian<double>& sq_hamiltonian, const size_t p, const size_t q) = 0;

    /**
     *  @param sq_hamiltonian       the current Hamiltonian
     *  @param p                    the index of spatial orbital 1
     *  @param q                    the index of spatial orbital 2
     *
     *  @return the angle for which the derivative of the scalar function after the Jacobi rotation is zero (and the second derivative is positive), using the current trigoniometric polynomial coefficients
     */
    virtual double calculateOptimalRotationAngle(const RSQHamiltonian<double>& sq_hamiltonian, const size_t p, const size_t q) const = 0;

    /**
     *  @param sq_hamiltonian               the current Hamiltonian
     *  @param jacobi_rot_par               the Jacobi rotation parameters
     * 
     *  @return the change in value for the scalar function if the given Jacobi rotation parameters would be used to rotate the given Hamiltonian
     */
    virtual double calculateScalarFunctionChange(const RSQHamiltonian<double>& sq_hamiltonian, const JacobiRotationParameters& jacobi_rot_par) const = 0;

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Jacobi-based orbital optimizer
     */
    virtual void prepareJacobiSpecificConvergenceChecking(const RSQHamiltonian<double>& sq_hamiltonian) = 0;


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param sq_hamiltonian           the current Hamiltonian
     * 
     *  @return a unitary matrix that will be used to rotate the current Hamiltonian into the next iteration
     */
    RTransformationMatrix<double> calculateNewRotationMatrix(const RSQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  @param sq_hamiltonian           the current Hamiltonian
     * 
     *  @return if the algorithm is considered to be converged
     */
    bool checkForConvergence(const RSQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
     */
    void prepareConvergenceChecking(const RSQHamiltonian<double>& sq_hamiltonian) override;


    // PUBLIC METHODS

    /**
     *  @param sq_hamiltonian           the Hamiltonian
     * 
     *  @return the optimal Jacobi rotation and the corresponding value for the scalar function that can be obtained when the Jacobi rotation would have taken place
     */
    std::pair<JacobiRotationParameters, double> calculateOptimalJacobiParameters(const RSQHamiltonian<double>& sq_hamiltonian);

    /**
     *  @return the comparer functor that is used to compare two pair_types
     */
    std::function<bool(const pair_type&, const pair_type&)> comparer() const;
};


}  // namespace GQCP
