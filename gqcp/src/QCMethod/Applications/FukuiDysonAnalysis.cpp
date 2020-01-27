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
#include "QCMethod/Applications/FukuiDysonAnalysis.hpp"

#include "Basis/transform.hpp"
#include "Processing/Properties/properties.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"



namespace GQCP {
namespace QCMethod {


/*
 * CONSTRUCTORS
 */

/**
 *  @param molecule                 the molecule that will be solved for
 *  @param basis_set                the basisset that should be used
 *  @param use_diis                 flags if one wants to use the DIIS for the RHF solve as opposed to Plain RHF
 */
FukuiDysonAnalysis::FukuiDysonAnalysis(const Molecule& molecule, const std::string& basis_set, const bool use_diis) : 
        molecule (molecule),
        spinor_basis (RSpinorBasis<double, GTOShell>(molecule, basis_set)),
        sq_hamiltonian (SQHamiltonian<double>::Molecular(this->spinor_basis, molecule)),  // in AO basis
        basis_set (basis_set)
{
    const auto K = this->spinor_basis.numberOfSpatialOrbitals();
    const auto N_P = this->molecule.numberOfElectrons()/2;

    // Define a molecule on which an RHF calculation is allowed
    auto restricted_molecule = this->molecule; 
    const bool open_shell_entry = (N_P*2 != this->molecule.numberOfElectrons());  // if the entry has an odd number of electrons, we create a molecule with an even number of electrons to allow RHF
    if (open_shell_entry) {
        restricted_molecule = Molecule(this->molecule.nuclearFramework(), this->molecule.totalNucleicCharge() - this->molecule.numberOfElectrons() +1);
    }

    // Perform DIIS or plain RHF given the flag in the constructor
    if (use_diis) {
        auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(restricted_molecule.numberOfElectrons(), this->sq_hamiltonian, this->spinor_basis.overlap().parameters());
        auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
        const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
        const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment).groundStateParameters();

        basisTransform(this->spinor_basis, this->sq_hamiltonian, rhf_parameters.coefficientMatrix());
    } else {
        auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), this->sq_hamiltonian, this->spinor_basis.overlap().parameters());
        auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
        const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
        const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

        basisTransform(this->spinor_basis, this->sq_hamiltonian, rhf_parameters.coefficientMatrix());
    }

    // In order to supply the correct arguments to the algorithm we choose different Fock spaces as fock_space1 should always have the highest occupation to fit the algorithm
    if (open_shell_entry) {
        this->fock_space1 = ProductFockSpace(K, N_P + 1, N_P);
        this->fock_space2 = ProductFockSpace(K, N_P, N_P);
    } else {
        this->fock_space1 = ProductFockSpace(K, N_P, N_P);
        this->fock_space2 = ProductFockSpace(K, N_P, N_P - 1); 
    }

    // Solve the FCI Hamiltonian eigenvalue problems with a Davidson algorithm
    FCI fci1 = FCI(fock_space1);
    FCI fci2 = FCI(fock_space2);
    
    CISolver ci_solver1 (fci1, sq_hamiltonian);
    CISolver ci_solver2 (fci2, sq_hamiltonian);

    VectorX<double> initial_g1 = fock_space1.HartreeFockExpansion();
    VectorX<double> initial_g2 = fock_space2.HartreeFockExpansion();
    DavidsonSolverOptions davidson_solver_options1 (initial_g1);
    DavidsonSolverOptions davidson_solver_options2 (initial_g2);
    ci_solver1.solve(davidson_solver_options1);
    ci_solver2.solve(davidson_solver_options2);

    // Retrieve the wavefunctions
    const auto wavefunction1 = ci_solver1.makeWavefunction();
    const auto wavefunction2 = ci_solver2.makeWavefunction();

    // Extract dyson coefficients
    this->dyson_coefficients = calculateDysonOrbitalCoefficients(wavefunction1, wavefunction2);
    
    // Calculate the Fukui matrix
    RDMCalculator rdm_calculator1 (fock_space1);
    rdm_calculator1.set_coefficients(wavefunction1.get_coefficients());
    RDMCalculator rdm_calculator2 (fock_space2);
    rdm_calculator2.set_coefficients(wavefunction2.get_coefficients());

    const auto onerdm1 = rdm_calculator1.calculate1RDMs().one_rdm;
    const auto onerdm2 = rdm_calculator2.calculate1RDMs().one_rdm;

    this->fukui_matrix = onerdm1 - onerdm2;

    // Diagonalize the Fukui matrix to retrieve Fukui naturals
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes1 (this->fukui_matrix);

    this->fukui_naturals = VectorX<double>(saes1.eigenvalues());
    this->fukui_vectors = OneRDM<double>(saes1.eigenvectors());
}


}  // namespace QCMethod
}  // namespace GQCP
