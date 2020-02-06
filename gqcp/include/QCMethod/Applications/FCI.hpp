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


#include "ONVBasis/ProductONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  A class that is a wrapper around solving the dense eigenvalue problem for the molecular Hamiltonian
 */
class FCI {
private:
    size_t N_alpha;  // the number of alpha electrons
    size_t N_beta;  // the number of beta electrons

    Molecule molecule;  // the molecule that will be solved for
    std::string basis_set;  // the basisset that should be used

    bool is_solved = false;
    bool use_davidson;
    double energy_solution;
    double s_squared;


public:
    // CONSTRUCTORS

    /**
     *  @param molecule             the molecule that will be solved for
     *  @param basis_set            the basisset that should be used
     *  @param num_alpha            the number of alpha electrons
     *  @param num_beta             the number of beta electrons
     */
    FCI(const Molecule& molecule, const std::string& basis_set, const size_t num_alpha, const size_t num_beta, const bool use_davidson);


    /**
     *  @param xyz_filename         the file that contains the molecule specification (coordinates in angstrom)
     *  @param basis_set            the basisset that should be used
     *  @param num_alpha            the number of alpha electrons
     *  @param num_beta             the number of beta electrons
     */
    FCI(const std::string& xyz_filename, const std::string& basis_set, const size_t num_alpha, const size_t num_beta, const bool use_davidson);


    // PUBLIC METHODS

    /**
     *  Solve the dense eigenvalue problem for the molecular Hamiltonian in the full ONV basis
     */
    void solve();

    /**
     *  @return the ground state FCI energy
     */
    double energy() const;

    /**
     *  @return the expectation value of the square of the spin angular momentum operator of the ground state
     */
    double spinSquared() const;
};


}  // namespace QCMethod
}  // namespace GQCP
