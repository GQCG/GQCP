// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef DOCINewtonOrbitalOptimizer_hpp
#define DOCINewtonOrbitalOptimizer_hpp

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "OrbitalOptimizationOptions.hpp"
#include "WaveFunction/WaveFunction.hpp"

#include <numopt.hpp>


namespace GQCP {


/**
 *  A class that performs gradient-and-Hessian-based orbital optimization for DOCI by sequentially
 *      - solving the DOCI eigenvalue problem
 *      - solving the Newton step to find the anti-Hermitian orbital rotation parameters
 *      - rotating the underlying spatial orbital basis
 */
class DOCINewtonOrbitalOptimizer {
private:
    GQCP::DOCI doci;  // the DOCI Hamiltonian builder
    GQCP::HamiltonianParameters ham_par;

    bool is_converged = false;
    std::vector<numopt::eigenproblem::Eigenpair> eigenpairs;  // eigenvalues and -vectors


public:
    // CONSTRUCTORS
    /**
     *  @param doci         the DOCI HamiltonianBuilder
     *  @param ham_par      the Hamiltonian parameters in an orthonormal basis
     */
    DOCINewtonOrbitalOptimizer(const GQCP::DOCI& doci, const GQCP::HamiltonianParameters& ham_par);


    // GETTERS
    const std::vector<numopt::eigenproblem::Eigenpair>& get_eigenpairs() const;
    const numopt::eigenproblem::Eigenpair& get_eigenpair(size_t index = 0) const;


    // PUBLIC METHODS
    /**
     *  Do the orbital optimization for DOCI
     *
     *  @param solver_options       solver options for the CI solver
     *  @param oo_options           options for the orbital optimization
     */
    void solve(numopt::eigenproblem::BaseSolverOptions& solver_options, const GQCP::OrbitalOptimizationOptions& oo_options=GQCP::OrbitalOptimizationOptions());

    /**
     *  @param index        the index of the index-th excited state
     *
     *  @return the index-th excited state after doing the OO-DOCI calculation
     */
    GQCP::WaveFunction get_wavefunction(size_t index = 0);
};


}  // namespace GQCP


#endif /* DOCINewtonOrbitalOptimizer_hpp */
