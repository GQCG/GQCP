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
#include "OrbitalOptimization/BaseOrbitalOptimizer.hpp"
#include "OrbitalOptimization/OrbitalRotationGenerators.hpp"
#include "Mathematical/Optimization/BaseHessianModifier.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Mathematical/Representation/SquareRankFourTensor.hpp"

#include <utility>


namespace GQCP {


/**
 *  An intermediate abstract class that should be derived from to implement a Newton-step based orbital optimization: the orbital gradient and Hessian are calculated through the DMs
 */
class NewtonOrbitalOptimizer : public BaseOrbitalOptimizer {
protected:
    std::shared_ptr<BaseHessianModifier> hessian_modifier;  // the modifier functor that should be used when an indefinite Hessian is encountered

    VectorX<double> gradient;
    SquareMatrix<double> hessian;


public:
    // CONSTRUCTORS

    /*
     *  @param hessian_modifier                 the modifier functor that should be used when an indefinite Hessian is encountered
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     */
    NewtonOrbitalOptimizer(std::shared_ptr<BaseHessianModifier> hessian_modifier, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128);


    // DESTRUCTOR
    virtual ~NewtonOrbitalOptimizer() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to calculate the first and second orbital derivatives, i.e. the orbital gradient and Hessian
     */
    virtual void prepareOrbitalDerivativesCalculation(const HamiltonianParameters<double>& ham_par) = 0;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital gradient as a matrix
     */
    virtual SquareMatrix<double> calculateGradientMatrix(const HamiltonianParameters<double>& ham_par) const = 0;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital Hessian as a tensor
     */
    virtual SquareRankFourTensor<double> calculateHessianTensor(const HamiltonianParameters<double>& ham_par) const = 0;

    /**
     *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the new full set orbital generators, including the redundant parameters
     */
    virtual OrbitalRotationGenerators calculateNewFullOrbitalGenerators(const HamiltonianParameters<double>& ham_par) const = 0;


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
     */
    virtual void prepareConvergenceChecking(const HamiltonianParameters<double>& ham_par) override;

    /**
     *  Determine if the algorithm has converged or not
     *  Specifically for the Newton-step based algorithms, this function
     *      - computes the gradient and checks its norm for convergence
     *      - if the gradient is zero, the Hessian is calculated and positive definiteness is checked
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return if the algorithm is considered to be converged
     */
    bool checkForConvergence(const HamiltonianParameters<double>& ham_par) const override;

    /**
     *  Produce a new rotation matrix by either
     *      - continuing in the direction of the i.e. the smallest (negative) eigenvalue
     *      - using the Newton step if it is well-defined
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return a unitary matrix that will be used to rotate the current Hamiltonian parameters into the next iteration
     */
    TransformationMatrix<double> calculateNewRotationMatrix(const HamiltonianParameters<double>& ham_par) const override;


    // PUBLIC METHODS

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital gradient as a vector. Matrix indices are converted to vector indices in the convention that p>q
     */
    VectorX<double> calculateGradientVector(const HamiltonianParameters<double>& ham_par) const;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital Hessian as a matrix
     */
    SquareMatrix<double> calculateHessianMatrix(const HamiltonianParameters<double>& ham_par) const;

    /**
     *  @return if a Newton step would be well-defined, i.e. the Hessian is positive definite
     */
    bool newtonStepIsWellDefined() const;

    /**
     *  If the Newton step is ill-defined, examine the Hessian and produce a new direction from it: the eigenvector that corresponds to the smallest (negative) eigenvalue of the Hessian
     * 
     *  @return the new direction from the Hessian if the Newton step is ill-defined
     */
    VectorX<double> directionFromIndefiniteHessian() const;

    /**
     *  Use gradient and Hessian information to determine a new direction for the 'free' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the new free orbital generators
     */
    OrbitalRotationGenerators calculateNewFreeOrbitalGenerators(const HamiltonianParameters<double>& ham_par)const;
};


}  // namespace GQCP
