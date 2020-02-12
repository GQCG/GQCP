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
#include "QCMethod/Applications/DOCIRHF.hpp"

#include "Basis/transform.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"


namespace GQCP {
namespace QCMethod {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param molecule             the molecule that will be solved for
 *  @param basis_set            the basisset that should be used
 *  @param use_davidson         indicate if one wants to use davidson to solve the eigenvalue problem (opposed to dense)
 */
DOCIRHF::DOCIRHF(const Molecule& molecule, const std::string& basis_set, const bool use_davidson) :
    molecule (molecule),
    basis_set (basis_set),
    use_davidson (use_davidson)
{}


/**
 *  @param xyz_filename         the file that contains the molecule specification (coordinates in angstrom)
 *  @param basis_set            the basisset that should be used
 *  @param use_davidson         indicate if one wants to use davidson to solve the eigenvalue problem (opposed to dense)
 */
DOCIRHF::DOCIRHF(const std::string& xyz_filename, const std::string& basis_set, const bool use_davidson) :
    DOCIRHF(Molecule::ReadXYZ(xyz_filename), basis_set, use_davidson)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the dense eigenvalue problem for the molecular Hamiltonian in the doubly occupied ONV basis
 */
void DOCIRHF::solve() {

    // Construct the molecular Hamiltonian in the RHF basis
    RSpinorBasis<double, GTOShell> spinor_basis (this->molecule, this->basis_set);
    const auto N_P = molecule.numberOfElectrons()/2;
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, this->molecule);  // in AO basis
    const size_t K = sq_hamiltonian.dimension();

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(this->molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();
    basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());

    // Set up DOCI
    SpinUnresolvedONVBasis fock_space (K, N_P);

    std::shared_ptr<BaseSolverOptions> solver_options;
    if (use_davidson) {
        solver_options = std::make_shared<DavidsonSolverOptions>(fock_space.hartreeFockExpansion());
    } else {
        solver_options = std::make_shared<DenseSolverOptions>();
    }

    DOCI doci (fock_space);

    CISolver solver (doci, sq_hamiltonian);
    solver.solve(*solver_options);
    double doci_energy = solver.get_eigenpair().get_eigenvalue();

    this->is_solved = true;
    double internuclear_repulsion_energy = Operator::NuclearRepulsion(this->molecule).value();
    this->energy_solution = doci_energy + internuclear_repulsion_energy;
    this->rhf_energy_solution = rhf_qc_structure.groundStateEnergy() + internuclear_repulsion_energy;
    this->T_total = spinor_basis.coefficientMatrix();
}


/**
 *  @return the DOCI energy
 */  
double DOCIRHF::energy() const {

    if (this->is_solved) {
        return this->energy_solution;
    } else {
        throw std::runtime_error("DOCIRHF::energy(): You are trying to get energy but the method hasn't been solved yet.");
    }
}


/**
 *  @return the RHF energy
 */  
double DOCIRHF::energy_rhf() const {

    if (this->is_solved) {
        return this->rhf_energy_solution;
    } else {
        throw std::runtime_error("DOCIRHF::energy_rhf(): You are trying to get energy but the method hasn't been solved yet.");
    }
}


/**
 *  @return the total transformation from atomic orbital basis to the RHF orbitals
 */
const TransformationMatrix<double>& DOCIRHF::transformationMatrix() const {

    if (this->is_solved) {
        return this->T_total;
    } else {
        throw std::runtime_error("DOCIRHF::transformationMatrix(): You are trying to get the total transformation matrix but the method hasn't been solved yet.");
    }
}


}  // namespace QCMethod
}  // namespace GQCP
