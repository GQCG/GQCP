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
#include "QCMethod/HF/GHF/GHFStabilityChecks.hpp"
#include "QCModel/HF/Stability/StabilityMatrices/GHFStabilityMatrix.hpp"

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

    auto environment = GQCP::GHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, S);

    // Perform a GHF SCF calculation
    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    auto ghf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // Initially, we don't know any of the stability properties.
    BOOST_CHECK(ghf_parameters.stabilityProperties().internal_stability == GQCP::Stability::unknown);
    BOOST_CHECK(ghf_parameters.stabilityProperties().external_stability == GQCP::Stability::unknown);

    // This method should be internally unstable, External should remain unknown.
    GQCP::QCMethod::GHFStabilityChecks<double>().internalStabilityCheck(ghf_parameters, sq_hamiltonian);

    BOOST_CHECK(ghf_parameters.stabilityProperties().internal_stability == GQCP::Stability::unstable);
    BOOST_CHECK(ghf_parameters.stabilityProperties().external_stability == GQCP::Stability::unknown);

    // This wavefunction should also be externally unstable. Internally, it should remain unstable.
    GQCP::QCMethod::GHFStabilityChecks<double>().externalStabilityCheck(ghf_parameters, sq_hamiltonian);

    BOOST_CHECK(ghf_parameters.stabilityProperties().internal_stability == GQCP::Stability::unstable);
    BOOST_CHECK(ghf_parameters.stabilityProperties().external_stability == GQCP::Stability::unstable);

    // Check that the stability properties can be printed
    ghf_parameters.stabilityProperties().print();
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
    GQCP::GTransformationMatrix<double> C_initial {6};
    // We take a randomly generated matrix (from @xdvriend's implementation) as the initial guess. This makes sure the off-diagonal spin-blocks are not zero blocks and helps the calculation to converge to a true GHF solution.
    // clang-format off
    C_initial << -0.3100721,  -0.15761163, -0.51612194, -0.38100148,  0.57090929, -0.37620802,
                 -0.00741269,  0.38801568, -0.25974834, -0.41043789, -0.67141074, -0.40332126,
                 -0.61961507,  0.18043708,  0.58367365,  0.17317687,  0.05464039, -0.45811451,
                  0.67031756,  0.28266352,  0.37079814, -0.23639173,  0.37758712, -0.3671939,
                  0.18059725, -0.8326703,   0.16282789, -0.03436191, -0.27832567, -0.41095738,
                  0.19477298,  0.13713633, -0.4018331,   0.77416187,  0.01572939, -0.42686445;
    // clang-format on
    GQCP::GHFSCFEnvironment<double> environment {N, sq_hamiltonian, S, C_initial};

    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 4000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    auto ghf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // Initially, we don't know any of the stability properties.
    BOOST_CHECK(ghf_parameters.stabilityProperties().internal_stability == GQCP::Stability::unknown);
    BOOST_CHECK(ghf_parameters.stabilityProperties().external_stability == GQCP::Stability::unknown);

    // This method should be internally unstable, External should remain unknown.
    GQCP::QCMethod::GHFStabilityChecks<double>().internalStabilityCheck(ghf_parameters, sq_hamiltonian);

    BOOST_CHECK(ghf_parameters.stabilityProperties().internal_stability == GQCP::Stability::stable);
    BOOST_CHECK(ghf_parameters.stabilityProperties().external_stability == GQCP::Stability::unknown);

    // This wavefunction should also be externally unstable. Internally, it should remain unstable.
    GQCP::QCMethod::GHFStabilityChecks<double>().externalStabilityCheck(ghf_parameters, sq_hamiltonian);

    BOOST_CHECK(ghf_parameters.stabilityProperties().internal_stability == GQCP::Stability::stable);
    BOOST_CHECK(ghf_parameters.stabilityProperties().external_stability == GQCP::Stability::stable);

    // Check that the stability properties can be printed
    ghf_parameters.stabilityProperties().print();
}

/**
 *  Starting from a core guess, the GHF SCF algorithm finds a solution that should be both internally and externally stable for the given system.
 *  This test checks whether the stability checks confirm this.
 * 
 *  The system of interest is a H2-chain, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H2_stability_test) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HChain(2, 1.0);  // H3-triangle, 1 bohr apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap().parameters();

    const auto sq_hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(g_spinor_basis, molecule);

    auto environment = GQCP::GHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, S);

    // Perform a GHF SCF calculation
    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    auto ghf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // Initially, we don't know any of the stability properties.
    BOOST_CHECK(ghf_parameters.stabilityProperties().internal_stability == GQCP::Stability::unknown);
    BOOST_CHECK(ghf_parameters.stabilityProperties().external_stability == GQCP::Stability::unknown);

    // This method should be internally unstable, External should remain unknown.
    GQCP::QCMethod::GHFStabilityChecks<double>().internalStabilityCheck(ghf_parameters, sq_hamiltonian);

    BOOST_CHECK(ghf_parameters.stabilityProperties().internal_stability == GQCP::Stability::stable);
    BOOST_CHECK(ghf_parameters.stabilityProperties().external_stability == GQCP::Stability::unknown);

    // This wavefunction should also be externally unstable. Internally, it should remain unstable.
    GQCP::QCMethod::GHFStabilityChecks<double>().externalStabilityCheck(ghf_parameters, sq_hamiltonian);

    BOOST_CHECK(ghf_parameters.stabilityProperties().internal_stability == GQCP::Stability::stable);
    BOOST_CHECK(ghf_parameters.stabilityProperties().external_stability == GQCP::Stability::stable);

    // Check that the stability properties can be printed
    ghf_parameters.stabilityProperties().print();
}
