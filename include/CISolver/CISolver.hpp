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
#ifndef GQCP_CISOLVER_HPP
#define GQCP_CISOLVER_HPP


#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "WaveFunction.hpp"

#include <numopt.hpp>


namespace GQCP {


/**
 *  Class which solves the CI eigenvalue problem and requires a HamiltonianBuilder and HamiltonianParameters
 *  so that it can find a set of eigenvalues and -vectors for the Hamiltonian (which cannot always be stored in memory)
 */
class CISolver {
private:
    HamiltonianBuilder* hamiltonian_builder;
    HamiltonianParameters hamiltonian_parameters;

    std::vector<numopt::eigenproblem::Eigenpair> eigenpairs;  // eigenvalues and -vectors

public:
    // CONSTRUCTOR
    /**
     *  Constructor given a @param hamiltonian_builder and @param hamiltonian_parameters
     */
    CISolver(HamiltonianBuilder& hamiltonian_builder, const HamiltonianParameters& hamiltonian_parameters);


    // GETTERS
    std::vector<numopt::eigenproblem::Eigenpair> get_eigenpairs() const { return this->eigenpairs; }
    numopt::eigenproblem::Eigenpair get_eigenpair(size_t index = 0) const { return this->eigenpairs[index]; }


    // PUBLIC METHODS
    /**
     *  solves the CI problem, setting the eigenpairs
     */
    void solve(numopt::eigenproblem::BaseSolverOptions& solver_options);

    /**
     *  @return WaveFunction instance after solving the CI problem for a given eigenvector at @param index
     */
    GQCP::WaveFunction get_wavefunction(size_t index = 0);
};


}  // namespace GQCP

#endif  // GQCP_CISOLVER_HPP
