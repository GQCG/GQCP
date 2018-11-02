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
#ifndef GQCP_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP
#define GQCP_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP



#include <memory>

#include "AOBasis.hpp"
#include "HamiltonianParameters.hpp"



namespace GQCP {


/**
 *  @param ao_basis     the AO basis in which the molecular Hamiltonian should be expressed
 *
 *  @return Hamiltonian parameters corresponding to the molecular Hamiltonian. The molecular Hamiltonian has
 *      - one-electron contributions:
 *          - kinetic
 *          - nuclear attraction
 *      - two-electron contributions:
 *          - Coulomb repulsion
 */
GQCP::HamiltonianParameters constructMolecularHamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis);


/**
 *  @param K        the number of orbitals
 *
 *  @return a set of random Hamiltonian parameters with values uniformly distributed between [-1,1]
 */
GQCP::HamiltonianParameters constructRandomHamiltonianParameters(size_t K);


/**
 *  @param fcidump_file     the name of the FCIDUMP file
 *
 *  @return Hamiltonian parameters corresponding to the contents of an FCIDUMP file
 */
GQCP::HamiltonianParameters readFCIDUMPFile(const std::string& fcidump_file);


/**
 *  @param upper_triagonal      an upper triagonal representation of the Hubbard hopping matrix
 *
 *  @return Hubbard Hamiltonian parameters generated from the Hubbard hopping matrix
 */
GQCP::HamiltonianParameters constructHubbardParameters(const Eigen::VectorXd& upper_triagonal);


}  // namespace GQCP


#endif  // GQCP_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP
