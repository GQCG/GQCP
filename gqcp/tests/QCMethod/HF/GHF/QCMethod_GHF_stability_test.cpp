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
    const auto S = g_spinor_basis.overlap().parameters();

    const auto sq_hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(g_spinor_basis, molecule);

    // Perform a GHF SCF calculation
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
    const auto S = g_spinor_basis.overlap().parameters();

    const auto sq_hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(g_spinor_basis, molecule);

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
    const auto S = g_spinor_basis.overlap().parameters();

    // Create a Hamiltonian in the AO basis.
    const auto sq_hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(g_spinor_basis, molecule);

    // Perform a GHF SCF calculation
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
