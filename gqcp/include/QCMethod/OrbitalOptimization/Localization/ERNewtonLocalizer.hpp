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


#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "QCMethod/OrbitalOptimization/NewtonOrbitalOptimizer.hpp"


namespace GQCP {


/**
 *  A class that localizes a set of orthonormal real spatial orbitals according to the maximization of the Edmiston-Ruedenberg localization index formulated as a minimization problem. The minimum is found using subsequent Newton steps.
 */
class ERNewtonLocalizer: public NewtonOrbitalOptimizer {
private:
    OrbitalSpace orbital_space;  // the occupied-virtual orbital space


public:
    // CONSTRUCTORS

    /**
     *  @param orbital_space                    the occupied-virtual orbital space
     *  @param hessian_modifier                 the modifier functor that should be used when an indefinite Hessian is encountered
     *  @param convergence_threshold            the threshold used to check for convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations that may be used to achieve convergence
     */
    ERNewtonLocalizer(const OrbitalSpace orbital_space, std::shared_ptr<BaseHessianModifier> hessian_modifier, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param sq_hamiltonian      the current Hamiltonian
     *
     *  @return the current orbital gradient of the Edmiston-Ruedenberg localization index as a matrix
     */
    SquareMatrix<double> calculateGradientMatrix(const RSQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  @param sq_hamiltonian      the current Hamiltonian
     *
     *  @return the current orbital Hessian of the Edmiston-Ruedenberg localization index as a tensor
     */
    SquareRankFourTensor<double> calculateHessianTensor(const RSQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
     * 
     *  @param sq_hamiltonian      the current Hamiltonian
     * 
     *  @return the new full set orbital generators, including the redundant parameters
     */
    OrbitalRotationGenerators calculateNewFullOrbitalGenerators(const RSQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
     */
    void prepareOrbitalDerivativesCalculation(const RSQHamiltonian<double>& sq_hamiltonian) override {}


    // PUBLIC METHODS

    /**
     *  @param sq_hamiltonian       the current Hamiltonian
     *  @param i                    the row of the gradient 'matrix'
     *  @param j                    the column of the gradient 'matrix'
     *
     *  @return the element (i,j) of the Edmiston-Ruedenberg localization index gradient
     */
    double calculateGradientMatrixElement(const RSQHamiltonian<double>& sq_hamiltonian, const size_t i, const size_t j) const;

    /**
     *  @param sq_hamiltonian       the current Hamiltonian
     *  @param i                    the first index of the Hessian 'tensor'
     *  @param j                    the second index of the Hessian 'tensor'
     *  @param k                    the third index of the Hessian 'tensor'
     *  @param l                    the fourth index of the Hessian 'tensor'
     *
     *  @return the element (i,j,k,l) of the Edmiston-Ruedenberg localization index Hessian
     */
    double calculateHessianTensorElement(const RSQHamiltonian<double>& sq_hamiltonian, const size_t i, const size_t j, const size_t k, const size_t l) const;
};


}  // namespace GQCP
