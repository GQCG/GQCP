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
#include "QCMethod/FCI.hpp"

namespace GQCP {
namespace QCMethod {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param xyz_filename         the file that contains the molecule specification (coordinates in angstrom)
 *  @param basis_set            the basisset that should be used
 *  @param num_alpha            the number of alpha electrons
 *  @param num_beta             the number of beta electrons
 */
FCI::FCI(const std::string xyz_filename, const std::string basis_set, const size_t num_alpha, const size_t num_beta) :
    xyz_filename (xyz_filename),
    basis_set (basis_set),
    N_alpha (num_alpha),
    N_beta (num_beta)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the dense eigenvalue problem for the molecular Hamiltonian in the full Fock space
 */
void FCI::solve() {

    // Construct the molecular Hamiltonian parameters
    auto molecule = Molecule::Readxyz(this->xyz_filename);
    auto mol_ham_par = HamiltonianParameters<double>::Molecular(molecule, this->basis_set);  // in the AO basis
    mol_ham_par.LowdinOrthonormalize();  // now in the Löwdin basis


    // Solve the FCI eigenvalue problem using the dense algorithm
    auto K = mol_ham_par.get_K();
    ProductFockSpace fock_space(K, this->N_alpha, this->N_beta);
    GQCP::FCI fci_builder (fock_space);

    CISolver fci_solver (fci_builder, mol_ham_par);
    DenseSolverOptions dense_solver_options;

    fci_solver.solve(dense_solver_options);
    this->is_solved = true;


    // Set the solution
    double fci_energy = fci_solver.get_eigenpair().get_eigenvalue();
    double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();

    this->energy_solution = fci_energy + internuclear_repulsion_energy;
}


/**
 *  @return the ground state FCI energy
 */
double FCI::energy() const {

    if (this->is_solved) {
        return this->energy_solution;
    } else {
        throw std::runtime_error("FCI::energy(): You are trying to get energy but the method hasn't been solved yet.");
    }
}


}  // namespace QCMethod
}  // namespace GQCP
