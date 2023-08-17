// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE "QCMethod_UHF_stability_test"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/UHF/UHF.hpp"
#include "QCMethod/HF/UHF/UHFSCFSolver.hpp"
#include "QCModel/HF/StabilityMatrices/UHFStabilityMatrices.hpp"


/**
 *  Starting from a core guess, the UHF SCF algorithm finds a solution that should be internally stable and externally unstable for the given system.
 *  This test checks whether the stability checks confirm this.
 *
 *  The system of interest is a H3-triangle, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_stability_test) {

    // Set up a spin orbital basis to obtain a second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.0);  // H3-triangle, 1 bohr apart.
    const auto N_alpha = molecule.numberOfElectronPairs() + (molecule.numberOfElectrons() - 2 * molecule.numberOfElectronPairs());
    const auto N_beta = molecule.numberOfElectronPairs();

    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto S = spinor_basis.overlap();

    const auto sq_hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In an AO basis.

    // Perform a UHF SCF calculation.
    auto environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, S);
    auto solver = GQCP::UHFSCFSolver<double>::Plain(1.0e-06, 3000);
    const auto qc_structure = GQCP::QCMethod::UHF<double>().optimize(solver, environment);
    auto uhf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need an unrestricted Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_unrestricted_mo = sq_hamiltonian.transformed(uhf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = uhf_parameters.calculateStabilityMatrices(hamiltonian_unrestricted_mo);

    // This method should be internally stable.
    const auto internal_stability = stability_matrices.isInternallyStable();
    BOOST_CHECK(internal_stability == true);

    // This wavefunction should be externally unstable.
    const auto external_stability = stability_matrices.isExternallyStable();
    BOOST_CHECK(external_stability == false);

    // We check both external instabilities separately as well.
    BOOST_CHECK(stability_matrices.isSpinUnconservedStable() == false);
    BOOST_CHECK(stability_matrices.isComplexConjugateStable() == true);

    // Check that the stability properties can be printed.
    stability_matrices.printStabilityDescription();
}

/**
 *  Starting from a core guess, the UHF SCF algorithm finds a solution that should be internally and externally unstable for the given system.
 *  This test checks whether the stability checks confirm this.
 *
 *  The system of interest is a H4-square, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H4_stability_test) {

    // Set up a spin orbital basis to obtain a second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(4, 1.0);  // H4-square, 1 bohr apart.
    const auto N_alpha = molecule.numberOfElectronPairs();
    const auto N_beta = molecule.numberOfElectronPairs();

    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "6-31G"};
    const auto S = spinor_basis.overlap();

    const auto sq_hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In an AO basis.

    // Perform a UHF SCF calculation.
    auto environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, S);
    auto solver = GQCP::UHFSCFSolver<double>::Plain(1.0e-06, 3000);
    const auto qc_structure = GQCP::QCMethod::UHF<double>().optimize(solver, environment);
    auto uhf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need an unrestricted Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_unrestricted_mo = sq_hamiltonian.transformed(uhf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = uhf_parameters.calculateStabilityMatrices(hamiltonian_unrestricted_mo);

    // This method should be internally unstable.
    const auto internal_stability = stability_matrices.isInternallyStable();
    BOOST_CHECK(internal_stability == false);

    // This wavefunction should be externally unstable.
    const auto external_stability = stability_matrices.isExternallyStable();
    BOOST_CHECK(external_stability == false);

    // We check both external instabilities separately as well.
    BOOST_CHECK(stability_matrices.isSpinUnconservedStable() == false);
    BOOST_CHECK(stability_matrices.isComplexConjugateStable() == false);

    // Check that the stability properties can be printed.
    stability_matrices.printStabilityDescription();
}


/**
 *  Starting from a core guess, the UHF SCF algorithm finds a solution that should be internally and externally unstable for the given system.
 *  This test checks whether the instability rotation matrix can rotate the guess towards the instability and find the lowest lying solution.
 *
 *  The system of interest is a H4-square, 2 Angstrom apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H4_stability_rotation) {

    // Set up a spin orbital basis to obtain a second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(4, 2 * 1.8897259886);  // H4-square, 1 bohr apart.
    const auto N_alpha = molecule.numberOfElectronPairs();
    const auto N_beta = molecule.numberOfElectronPairs();

    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "6-31G"};
    const auto S = spinor_basis.overlap();

    const auto sq_hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In an AO basis.

    // Perform a UHF SCF calculation.
    auto environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, S);
    auto solver = GQCP::UHFSCFSolver<double>::Plain(1.0e-06, 3000);
    const auto qc_structure = GQCP::QCMethod::UHF<double>().optimize(solver, environment);
    auto uhf_parameters = qc_structure.groundStateParameters();
    auto nuc_rep = GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();

    // We can now check the stability of the ground state parameters.
    // For this we need an unrestricted Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_unrestricted_mo = sq_hamiltonian.transformed(uhf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = uhf_parameters.calculateStabilityMatrices(hamiltonian_unrestricted_mo);

    const auto rotation = stability_matrices.instabilityRotationMatrix(N_alpha, N_beta, 6, 6);
    const auto new_guess = uhf_parameters.expansion().transformed(rotation);

    // Perform a new GHF calculation. Label `_rotated` is used to denote the variables AFTER rotation towards steepest descent.
    GQCP::UHFSCFEnvironment<double> environment_rotated {N_alpha, N_beta, sq_hamiltonian, S, new_guess};
    auto solver_rotated = GQCP::UHFSCFSolver<double>::Plain(1.0e-06, 3000);
    const auto qc_structure_rotated = GQCP::QCMethod::UHF<double>().optimize(solver_rotated, environment_rotated);
    auto uhf_parameters_rotated = qc_structure_rotated.groundStateParameters();

    // Check that the new energy corresponds with the stable solution.
    BOOST_CHECK((qc_structure_rotated.groundStateEnergy() + nuc_rep) - (-2.0066898106390494) < 1e-6);

    // Check the stability of the new GHF parameters.
    // We can now check the stability of the new ground state parameters.
    // For this we need a generalized Hamiltonian in the orthonormal MO basis of the new parameters.
    const auto hamiltonian_generalized_rotated = sq_hamiltonian.transformed(uhf_parameters_rotated.expansion());

    // Calculate the new stability matrices.
    const auto stability_matrices_rotated = uhf_parameters_rotated.calculateStabilityMatrices(hamiltonian_generalized_rotated);

    // Check that the stability properties can be printed.
    stability_matrices_rotated.printStabilityDescription();
}


/**
 *  Starting from a core guess, the complex valued UHF SCF algorithm finds a solution that should be internally and externally unstable for the given system.
 *  This test checks whether the instability rotation matrix can rotate the guess towards the instability and find the lowest lying solution.
 *
 *  The system of interest is a H4-square, 2 Angstrom apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H4_stability_rotation_complex) {

    // Set up a spin orbital basis to obtain a second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(4, 2 * 1.8897259886);  // H4-square, 1 bohr apart.
    const auto N_alpha = molecule.numberOfElectronPairs();
    const auto N_beta = molecule.numberOfElectronPairs();

    const GQCP::USpinOrbitalBasis<GQCP::complex, GQCP::GTOShell> spinor_basis {molecule, "6-31G"};
    const auto S = spinor_basis.overlap();

    const auto sq_hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In an AO basis.

    // Perform a UHF SCF calculation.
    auto environment = GQCP::UHFSCFEnvironment<GQCP::complex>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, S);
    auto solver = GQCP::UHFSCFSolver<GQCP::complex>::Plain(1.0e-06, 3000);
    const auto qc_structure = GQCP::QCMethod::UHF<GQCP::complex>().optimize(solver, environment);
    auto uhf_parameters = qc_structure.groundStateParameters();
    auto nuc_rep = GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();

    // We can now check the stability of the ground state parameters.
    // For this we need an unrestricted Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_unrestricted_mo = sq_hamiltonian.transformed(uhf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = uhf_parameters.calculateStabilityMatrices(hamiltonian_unrestricted_mo);

    const auto rotation = stability_matrices.instabilityRotationMatrix(N_alpha, N_beta, 6, 6);
    const auto new_guess = uhf_parameters.expansion().transformed(rotation);

    // Perform a new GHF calculation. Label `_rotated` is used to denote the variables AFTER rotation towards steepest descent.
    GQCP::UHFSCFEnvironment<GQCP::complex> environment_rotated {N_alpha, N_beta, sq_hamiltonian, S, new_guess};
    auto solver_rotated = GQCP::UHFSCFSolver<GQCP::complex>::Plain(1.0e-06, 3000);
    const auto qc_structure_rotated = GQCP::QCMethod::UHF<GQCP::complex>().optimize(solver_rotated, environment_rotated);
    auto uhf_parameters_rotated = qc_structure_rotated.groundStateParameters();

    // Check that the new energy corresponds with the stable solution.
    BOOST_CHECK((qc_structure_rotated.groundStateEnergy().real() + nuc_rep) - (-2.0066898106390494) < 1e-6);

    // Check the stability of the new GHF parameters.
    // We can now check the stability of the new ground state parameters.
    // For this we need a generalized Hamiltonian in the orthonormal MO basis of the new parameters.
    const auto hamiltonian_generalized_rotated = sq_hamiltonian.transformed(uhf_parameters_rotated.expansion());

    // Calculate the new stability matrices.
    const auto stability_matrices_rotated = uhf_parameters_rotated.calculateStabilityMatrices(hamiltonian_generalized_rotated);

    // Check that the stability properties can be printed.
    stability_matrices_rotated.printStabilityDescription();
}

/**
 *  Starting from a core guess, the UHF SCF algorithm finds a solution that should be internally unstable.
 *  This test checks whether the stability checks confirm this. It also rotates (twice) towards the lowest eigenvector to finally find the stable solution.
 *
 *  The system of interest is a H3-triangle, 1 bohr apart and the reference implementation was done in PySCF.
 */
BOOST_AUTO_TEST_CASE(LiF_internal_instability) {

    // Set up a spin orbital basis to obtain a second-quantized molecular Hamiltonian.
    const auto left = GQCP::Nucleus(3, 0, 0, 0);
    const auto right = GQCP::Nucleus(9, 0, 0, 2 * 1.88973);
    const auto molecule = GQCP::Molecule({left, right}, 0);  // H3-triangle, 1 bohr apart.
    const auto N_alpha = molecule.numberOfElectronPairs();
    const auto N_beta = molecule.numberOfElectronPairs();
    const auto nuc_rep = GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();

    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto S = spinor_basis.overlap();

    const auto sq_hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In an AO basis.

    // Perform a UHF SCF calculation.
    auto environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, S);
    auto solver = GQCP::UHFSCFSolver<double>::Plain(1.0e-06, 3000);
    const auto qc_structure = GQCP::QCMethod::UHF<double>().optimize(solver, environment);
    auto uhf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need an unrestricted Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_unrestricted_mo = sq_hamiltonian.transformed(uhf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = uhf_parameters.calculateStabilityMatrices(hamiltonian_unrestricted_mo);

    // This method should be internally stable.
    stability_matrices.printStabilityDescription();

    // Check energy and internal stability.
    BOOST_CHECK((qc_structure.groundStateEnergy() + nuc_rep) - (-105.295362733027) < 1e-6);
    BOOST_CHECK(stability_matrices.isInternallyStable() == false);

    // follow instability.
    const auto K_alpha = spinor_basis.numberOfSpinors() / 2;
    const auto K_beta = K_alpha;
    const auto Va = K_alpha - N_alpha;
    const auto Vb = K_beta - N_beta;

    const auto rotation_matrix = stability_matrices.instabilityRotationMatrix(N_alpha, N_beta, Va, Vb);
    const auto new_guess = uhf_parameters.expansion().rotated(rotation_matrix);

    auto environment_2 = GQCP::UHFSCFEnvironment<double>(N_alpha, N_beta, sq_hamiltonian, S, new_guess);
    const auto qc_structure_2 = GQCP::QCMethod::UHF<double>().optimize(solver, environment_2);

    auto uhf_parameters_2 = qc_structure_2.groundStateParameters();
    const auto hamiltonian_unrestricted_mo_2 = sq_hamiltonian.transformed(uhf_parameters_2.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices_2 = uhf_parameters_2.calculateStabilityMatrices(hamiltonian_unrestricted_mo_2);

    // This method should be internally stable.
    stability_matrices_2.printStabilityDescription();

    // Check energy and internal stability.
    BOOST_CHECK((qc_structure_2.groundStateEnergy() + nuc_rep) - (-105.317595983304) < 1e-6);
    BOOST_CHECK(stability_matrices_2.isInternallyStable() == false);

    // follow instability.
    const auto rotation_matrix_2 = stability_matrices_2.instabilityRotationMatrix(N_alpha, N_beta, Va, Vb);
    const auto new_guess_2 = uhf_parameters_2.expansion().rotated(rotation_matrix_2);

    auto environment_3 = GQCP::UHFSCFEnvironment<double>(N_alpha, N_beta, sq_hamiltonian, S, new_guess_2);
    const auto qc_structure_3 = GQCP::QCMethod::UHF<double>().optimize(solver, environment_3);

    auto uhf_parameters_3 = qc_structure_3.groundStateParameters();
    const auto hamiltonian_unrestricted_mo_3 = sq_hamiltonian.transformed(uhf_parameters_3.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices_3 = uhf_parameters_3.calculateStabilityMatrices(hamiltonian_unrestricted_mo_3);

    // This method should be internally stable.
    stability_matrices_3.printStabilityDescription();

    // Check energy and internal stability.
    BOOST_CHECK((qc_structure_3.groundStateEnergy() + nuc_rep) - (-105.325920420114) < 1e-6);
    BOOST_CHECK(stability_matrices_3.isInternallyStable() == true);
}
