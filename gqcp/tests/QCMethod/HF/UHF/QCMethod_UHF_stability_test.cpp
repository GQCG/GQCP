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

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto S = spinor_basis.overlap().parameters();

    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.

    // Perform a UHF SCF calculation.
    auto environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, S);
    auto solver = GQCP::UHFSCFSolver<double>::Plain(1.0e-06, 3000);
    const auto qc_structure = GQCP::QCMethod::UHF<double>().optimize(solver, environment);
    auto uhf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need an unrestricted Hamiltonian in the orthonormal MO basis.
    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> unrestricted_basis {molecule, "STO-3G"};
    const auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(unrestricted_basis, molecule);
    const auto hamiltonian_unrestricted = usq_hamiltonian.transformed(uhf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = uhf_parameters.calculateStabilityMatrices(hamiltonian_unrestricted);

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

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "6-31G"};
    const auto S = spinor_basis.overlap().parameters();

    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.

    // Perform a UHF SCF calculation.
    auto environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, S);
    auto solver = GQCP::UHFSCFSolver<double>::Plain(1.0e-06, 3000);
    const auto qc_structure = GQCP::QCMethod::UHF<double>().optimize(solver, environment);
    auto uhf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need an unrestricted Hamiltonian in the orthonormal MO basis.
    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> unrestricted_basis {molecule, "6-31G"};
    const auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(unrestricted_basis, molecule);
    const auto hamiltonian_unrestricted = usq_hamiltonian.transformed(uhf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = uhf_parameters.calculateStabilityMatrices(hamiltonian_unrestricted);

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
