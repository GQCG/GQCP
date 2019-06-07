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
#ifndef BASEORBITALOPTIMIZER_HPP
#define BASEORBITALOPTIMIZER_HPP


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "OrbitalOptimization/OrbitalOptimizationOptions.hpp"


namespace GQCP {


/**
 *  A base class for orbital optimizers. Due to the generality of the nature of the orbital optimization problem, the main algorithm (see solve()) is implemented inside this base class
 */
class BaseOrbitalOptimizer {
protected:
    OrbitalOptimizationOptions oo_options;


public:
    // CONSTRUCTORS
    /**
     *  @param oo_options       the orbital optimization options that should be used for the orbital optimization algorithm
     */
    BaseOrbitalOptimizer(const OrbitalOptimizationOptions& oo_options);


    // PUBLIC PURE VIRTUAL METHODS
    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return if the algorithm is considered to be converged
     */
    virtual bool checkForConvergence(const HamiltonianParameters<double>& ham_par) = 0;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return a unitary matrix that will be used to rotate the current Hamiltonian parameters into the next iteration
     */
    virtual SquareMatrix<double> calculateNewRotationMatrix(const HamiltonianParameters<double>& ham_par) = 0;


    // PUBLIC METHODS
    /**
     *  Optimize the Hamiltonian parameters by subsequently
     *      - checking for convergence (see checkForConvergence())
     *      - rotating the Hamiltonian parameters with a newly found rotation matrix (see calculateNewRotationMatrix())
     * 
     *  @param ham_par      the initial (guess for the) Hamiltonian parameters
     */
    void optimize(HamiltonianParameters<double>& ham_par);
};


}  // namespace GQCP



#endif  /* BASEORBITALOPTIMIZER_HPP */
