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
#ifndef GQCP_QCMETHOD_DOCINEWTONORBITALOPTIMIZER_HPP
#define GQCP_QCMETHOD_DOCINEWTONORBITALOPTIMIZER_HPP


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "RDM/RDMCalculator.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  A class that is a wrapper around solving the dense eigenvalue problem for the molecular Hamiltonian
 */
class DOCINewtonOrbitalOptimizer {
private:
    size_t N_P;

    std::string xyz_filename;  // the file that contains the molecule specification (coordinates in angstrom)
    std::string basis_set;  // the basisset that should be used

    bool is_solved = false;
    double energy_solution;


public:
    // CONSTRUCTORS

    /**
     *  @param xyz_filename         the file that contains the molecule specification (coordinates in angstrom)
     *  @param basis_set            the basisset that should be used
     *  @param num_alpha            the number of alpha electrons
     *  @param num_beta             the number of beta electrons
     */
    DOCINewtonOrbitalOptimizer(const std::string xyz_filename, const std::string basis_set, const bool use_davidson, const bool localize);


    // PUBLIC METHODS

    /**
     *  Solve the dense eigenvalue problem for the molecular Hamiltonian in the full Fock space
     */
    void solve();

    /**
     *  @return the ground state FCI energy
     */
    double energy() const;
};


}  // namespace QCMethod
}  // namespace GQCP


#endif  // GQCP_QCMETHOD_DOCINEWTONORBITALOPTIMIZER_HPP
