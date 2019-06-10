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
#ifndef GQCP_JACOBIORBITALOPTIMIZER_HPP
#define GQCP_JACOBIORBITALOPTIMIZER_HPP


#include "OrbitalOptimization/BaseOrbitalOptimizer.hpp"
#include "OrbitalOptimization/JacobiRotationParameters.hpp"

#include <utility>


namespace GQCP {


class JacobiOrbitalOptimizer : public BaseOrbitalOptimizer {
protected:
    size_t dim; // the dimension of the orbital space that should be scanned. The valid orbital indices then are 0 ... dim (not included)

    std::pair<JacobiRotationParameters, double> optimal_jacobi_with_scalar;  // holds the optimal Jacobi parameters and the corresponding value for the scalar function trying to optimize


public:
    // CONSTRUCTORS
    
    /**
     *  @param dim              the dimension of the orbital space that should be scanned. The valid orbital indices then are 0 ... dim (not included)
     *  @param oo_options       the options for orbital optimization
     */
    JacobiOrbitalOptimizer(const size_t dim, const OrbitalOptimizationOptions& oo_options);


    // PUBLIC PURE VIRTUAL METHODS
    /**
     *  @param ham_par      the Hamiltonian parameters
     * 
     *  @return the value of the scalar function that should be optimized
     */
    virtual double calculateScalarFunction(const HamiltonianParameters<double>& ham_par) = 0;

    /**
     *  Calculate the trigoniometric polynomial coefficients for the given Jacobi rotation
     *
     *  @param p            the index of spatial orbital p
     *  @param q            the index of spatial orbital q
     */
    virtual void calculateJacobiCoefficients(const HamiltonianParameters<double>& ham_par, const size_t p, const size_t q) = 0;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     *  @param p            the index of spatial orbital 1
     *  @param q            the index of spatial orbital 2
     *
     *  @return the angle for which the derivative of the scalar function after the Jacobi rotation is zero (and the second derivative is positive), using the current trigoniometric polynomial coefficients
     */
    virtual double calculateOptimalRotationAngle(const HamiltonianParameters<double>& ham_par, const size_t p, const size_t q) = 0;

    /**
     *  @param ham_par              the current Hamiltonian parameters
     *  @param jacobi_rot_par       the Jacobi rotation parameters
     * 
     *  @return the value of the scalar function if the given Jacobi rotation parameters would be used to rotate the given Hamiltonian parameters
     */
    virtual double calculateScalarFunctionAfterJacobiRotation(const HamiltonianParameters<double>& ham_par, const JacobiRotationParameters& jacobi_rot_par) = 0;


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return if the algorithm is considered to be converged
     */
    bool checkForConvergence(const HamiltonianParameters<double>& ham_par) override;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return a unitary matrix that will be used to rotate the current Hamiltonian parameters into the next iteration
     */
    SquareMatrix<double> calculateNewRotationMatrix(const HamiltonianParameters<double>& ham_par) override;


    // PUBLIC METHODS

    /**
     *  @param ham_par      the Hamiltonian parameters
     * 
     *  @return the optimal Jacobi rotation parameters and the corresponding value for the scalar function that can be obtained when the Jacobi rotation would have taken place
     */
    std::pair<JacobiRotationParameters, double> calculateOptimalJacobiParameters(const HamiltonianParameters<double>& ham_par);
};



}  // GQCP


#endif  // GQCP_JACOBIORBITALOPTIMIZER_HPP