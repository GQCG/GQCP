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

#define BOOST_TEST_MODULE "QCMethod_GHF_stability_test"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/GHF/GHF.hpp"
#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"
#include "QCModel/HF/StabilityMatrices/GHFStabilityMatrices.hpp"

/**
 *  Starting from a core guess, the GHF SCF algorithm finds a solution that should be both internally and externally unstable for the given system.
 *  This test checks whether the stability checks confirm this.
 *
 *  The system of interest is a H3-triangle, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_stability_test_1) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.0);  // H3-triangle, 1 bohr apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Perform a GHF SCF calculation.
    auto environment = GQCP::GHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, S);
    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    auto ghf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need a generalized Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_generalized = sq_hamiltonian.transformed(ghf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = ghf_parameters.calculateStabilityMatrices(hamiltonian_generalized);


    // This method should be internally unstable.
    const auto internal_stability = stability_matrices.isInternallyStable();
    BOOST_CHECK(internal_stability == false);

    // This wavefunction should also be externally unstable.
    const auto external_stability = stability_matrices.isExternallyStable();
    BOOST_CHECK(external_stability == false);

    // Check that the stability properties can be printed
    stability_matrices.printStabilityDescription();
}


/**
 *  Starting from a random guess, the GHF SCF algorithm finds a solution that should be both internally and externally stable for the given system.
 *  This test checks whether the stability checks confirm this.
 *
 *  The system of interest is a H3-triangle, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_stability_test_2) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.0);  // H3-triangle, 1 bohr apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Create a solver and associated environment and let the QCMethod do its job.
    GQCP::SquareMatrix<double> C_initial_matrix {6};
    // We take a randomly generated matrix (from @xdvriend's implementation) as the initial guess. This makes sure the off-diagonal spin-blocks are not zero blocks and helps the calculation to converge to a true GHF solution.
    // clang-format off
    C_initial_matrix << -0.3100721,  -0.15761163, -0.51612194, -0.38100148,  0.57090929, -0.37620802,
                        -0.00741269,  0.38801568, -0.25974834, -0.41043789, -0.67141074, -0.40332126,
                        -0.61961507,  0.18043708,  0.58367365,  0.17317687,  0.05464039, -0.45811451,
                         0.67031756,  0.28266352,  0.37079814, -0.23639173,  0.37758712, -0.3671939,
                         0.18059725, -0.8326703,   0.16282789, -0.03436191, -0.27832567, -0.41095738,
                         0.19477298,  0.13713633, -0.4018331,   0.77416187,  0.01572939, -0.42686445;
    // clang-format on
    GQCP::GTransformation<double> C_initial {C_initial_matrix};

    GQCP::GHFSCFEnvironment<double> environment {N, sq_hamiltonian, S, C_initial};
    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 4000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    auto ghf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need a generalized Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_generalized = sq_hamiltonian.transformed(ghf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = ghf_parameters.calculateStabilityMatrices(hamiltonian_generalized);


    // This method should be internally stable.
    const auto internal_stability = stability_matrices.isInternallyStable();
    BOOST_CHECK(internal_stability == true);

    // This wavefunction should also be externally stable.
    const auto external_stability = stability_matrices.isExternallyStable();
    BOOST_CHECK(external_stability == true);

    // Check that the stability properties can be printed
    stability_matrices.printStabilityDescription();
}

/**
 *  Starting from a core guess, the GHF SCF algorithm finds a solution that should be both internally and externally stable for the given system.
 *  This test checks whether the stability checks confirm this.
 *
 *  The system of interest is a H2-chain, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_stability_test_3) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.8897);  // H3-triangle, 1 angstrom apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "6-31G"};
    const auto S = g_spinor_basis.overlap();

    // Create a Hamiltonian in the AO basis.
    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Perform a GHF SCF calculation.
    auto environment = GQCP::GHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, S);
    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-05, 5000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    auto ghf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need a generalized Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_generalized = sq_hamiltonian.transformed(ghf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = ghf_parameters.calculateStabilityMatrices(hamiltonian_generalized);

    // This method should be internally stable.
    const auto internal_stability = stability_matrices.isInternallyStable();
    BOOST_CHECK(internal_stability == false);

    // This wavefunction should also be externally stable.
    const auto external_stability = stability_matrices.isExternallyStable();
    BOOST_CHECK(external_stability == false);

    // Check that the stability properties can be printed
    stability_matrices.printStabilityDescription();
}

/**
 *  This test checks whether the internal instability found in the real GHF calculation for H3 can be followed into a stable solution for that system.
 *
 *  The system of interest is a H3-triangle, 1 Angstrom apart and the reference implementation was done by @xdvriend in the "Following GHF instabilities" example on the GQCP website (https://gqcg.github.io/GQCP/examples/Following-internal-GHF-instabilities.html).
 */
BOOST_AUTO_TEST_CASE(H3_stability_rotation) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.889);  // H3-triangle, 1 Angstrom apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Perform a GHF SCF calculation.
    auto environment = GQCP::GHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, S);
    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    auto ghf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need a generalized Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_generalized = sq_hamiltonian.transformed(ghf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = ghf_parameters.calculateStabilityMatrices(hamiltonian_generalized);

    const auto occ = ghf_parameters.numberOfElectrons();
    const auto virt = ghf_parameters.numberOfSpinors() - occ;
    const auto rotation = stability_matrices.instabilityRotationMatrix(occ, virt);

    // Set up a reference rotation matrix, taken from the python example on the website.
    GQCP::SquareMatrix<double> reference {6};
    // clang-format off
    reference <<  0.984538137,  0.0,        0.0,         0.0,         0.0,         0.175170363,
                  0.0,          1.00000000, 0.0,         0.0,         0.0,         0.0,
                  0.0,          0.0,        0.553382907, 0.0,        -0.832926982, 0.0,
                  0.0,          0.0,        0.0,         1.00000000,  0.0,         0.0, 
                  0.0,          0.0,        0.832926982, 0.0,         0.553382907, 0.0,
                 -0.175170363,  0.0,        0.0,         0.0,         0.0,         0.984538137;
    // clang-format on

    // Check the calculated matrix versus the reference.
    BOOST_CHECK(rotation.matrix().isApprox(reference, 1.0e-06));

    // Transform the guess with the newly found rotation matrix.
    const auto new_guess = ghf_parameters.expansion().transformed(rotation);

    // Perform a new GHF calculation. Label `_rotated` is used to denote the variables AFTER rotation towards steepest descent.
    GQCP::GHFSCFEnvironment<double> environment_rotated {N, sq_hamiltonian, S, new_guess};
    auto solver_rotated = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 4000);
    const auto qc_structure_rotated = GQCP::QCMethod::GHF<double>().optimize(solver_rotated, environment_rotated);
    auto ghf_parameters_rotated = qc_structure_rotated.groundStateParameters();

    // Check the stability of the new GHF parameters.
    // We can now check the stability of the new ground state parameters.
    // For this we need a generalized Hamiltonian in the orthonormal MO basis of the new parameters.
    const auto hamiltonian_generalized_rotated = sq_hamiltonian.transformed(ghf_parameters_rotated.expansion());

    // Calculate the new stability matrices.
    const auto stability_matrices_rotated = ghf_parameters_rotated.calculateStabilityMatrices(hamiltonian_generalized_rotated);

    // This method should now be internally stable.
    const auto internal_stability = stability_matrices_rotated.isInternallyStable();
    BOOST_CHECK(internal_stability == true);

    // This wavefunction should now also be externally stable.
    const auto external_stability = stability_matrices_rotated.isExternallyStable();
    BOOST_CHECK(external_stability == true);

    // Check that the stability properties can be printed & match the stability conditions tested previously.
    stability_matrices_rotated.printStabilityDescription();

    // Check that the electronic energy now matches that of the stable GHF solution.
    BOOST_CHECK(qc_structure_rotated.groundStateEnergy() - (-2.9284445024360175) < 1e-8);
}


/**
 *  This test checks whether the a correct rotation matrix for complex valued GHF calculations for H3 can be followed into a stable solution for that system.
 *
 *  The system of interest is a H3-triangle, 1 Angstrom apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_stability_rotation_complex) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.889);  // H3-triangle, 1 Angstrom apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Perform a GHF SCF calculation.
    auto environment = GQCP::GHFSCFEnvironment<GQCP::complex>::WithCoreGuess(N, sq_hamiltonian, S);
    auto solver = GQCP::GHFSCFSolver<GQCP::complex>::Plain(1.0e-08, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<GQCP::complex>().optimize(solver, environment);
    auto ghf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need a generalized Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_generalized = sq_hamiltonian.transformed(ghf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = ghf_parameters.calculateStabilityMatrices(hamiltonian_generalized);

    const auto occ = ghf_parameters.numberOfElectrons();
    const auto virt = ghf_parameters.numberOfSpinors() - occ;
    const auto rotation = stability_matrices.instabilityRotationMatrix(occ, virt);

    // Transform the guess with the newly found rotation matrix.
    const auto new_guess = ghf_parameters.expansion().transformed(rotation);

    // Perform a new GHF calculation. Label `_rotated` is used to denote the variables AFTER rotation towards steepest descent.
    GQCP::GHFSCFEnvironment<GQCP::complex> environment_rotated {N, sq_hamiltonian, S, new_guess};
    auto solver_rotated = GQCP::GHFSCFSolver<GQCP::complex>::Plain(1.0e-08, 4000);
    const auto qc_structure_rotated = GQCP::QCMethod::GHF<GQCP::complex>().optimize(solver_rotated, environment_rotated);
    auto ghf_parameters_rotated = qc_structure_rotated.groundStateParameters();

    // Check the stability of the new GHF parameters.
    // We can now check the stability of the new ground state parameters.
    // For this we need a generalized Hamiltonian in the orthonormal MO basis of the new parameters.
    const auto hamiltonian_generalized_rotated = sq_hamiltonian.transformed(ghf_parameters_rotated.expansion());

    // Calculate the new stability matrices.
    const auto stability_matrices_rotated = ghf_parameters_rotated.calculateStabilityMatrices(hamiltonian_generalized_rotated);

    // This method should now be internally stable.
    const auto internal_stability = stability_matrices_rotated.isInternallyStable();
    BOOST_CHECK(internal_stability == true);

    // Check that the stability properties can be printed & match the stability conditions tested previously.
    stability_matrices_rotated.printStabilityDescription();

    // Check that the electronic energy now matches that of the stable GHF solution.
    BOOST_CHECK(qc_structure_rotated.groundStateEnergy().real() - (-2.9284445024360175) < 1e-8);
}
