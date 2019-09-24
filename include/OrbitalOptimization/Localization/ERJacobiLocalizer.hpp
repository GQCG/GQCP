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


#include "OrbitalOptimization/JacobiOrbitalOptimizer.hpp"


namespace GQCP {


/**
 *  A class that localizes a set of orthonormal orbitals according to the maximization of the Edmiston-Ruedenberg localization index formulated as a minimization problem. The minimum is found using subsequent Jacobi rotations
 */
class ERJacobiLocalizer : public JacobiOrbitalOptimizer {
private:
    double A=0.0, B=0.0, C=0.0;  // the Jacobi rotation coefficients


public:
    // CONSTRUCTORS
    using JacobiOrbitalOptimizer::JacobiOrbitalOptimizer;  // inherit base constructor


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence
     */
    void prepareJacobiSpecificConvergenceChecking(const SQHamiltonian<double>& sq_hamiltonian) override {}

    /**
     *  Calculate the trigoniometric polynomial coefficients for the given Jacobi rotation indices
     *
     *  @param i            the index of spatial orbital 1
     *  @param j            the index of spatial orbital 2
     */
    void calculateJacobiCoefficients(const SQHamiltonian<double>& sq_hamiltonian, const size_t i, const size_t j) override;

    /**
     *  @param sq_hamiltonian       the current Hamiltonian
     *  @param i                    the index of spatial orbital 1
     *  @param j                    the index of spatial orbital 2
     *
     *  @return the angle for which the derivative of the scalar function after the Jacobi rotation is zero (and the second derivative is positive), using the current trigoniometric polynomial coefficients
     */
    double calculateOptimalRotationAngle(const SQHamiltonian<double>& sq_hamiltonian, const size_t i, const size_t j) const override;

    /**
     *  @param sq_hamiltonian           the current Hamiltonian
     *  @param jacobi_rot_par           the Jacobi rotation parameters
     * 
     *  @return the change in the value of the scalar function (i.e. minus the ER localization index) if the given Jacobi rotation parameters would be used to rotate the given Hamiltonian
     */
    double calculateScalarFunctionChange(const SQHamiltonian<double>& sq_hamiltonian, const JacobiRotationParameters& jacobi_rot_par) const override;
};


}  // namespace GQCP
