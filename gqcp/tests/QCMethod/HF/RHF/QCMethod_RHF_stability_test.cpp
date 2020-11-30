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

#define BOOST_TEST_MODULE "QCMethod_RHF_stability_test"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCModel/HF/StabilityMatrices/RHFStabilityMatrices.hpp"


/**
 *  The RHF solution of H_2O should be stable both internally and externally.
 *  This test checks whether the RHF stability analysis confirms this.
 */
BOOST_AUTO_TEST_CASE(h2o_sto3g_stability) {

    // Do our own RHF calculation
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {water, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, water);  // in an AO basis
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(water.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const auto qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    auto rhf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need a restricted Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_restricted = sq_hamiltonian.transformed(rhf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = rhf_parameters.calculateStabilityMatrices(hamiltonian_restricted);

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
 *  The RHF solution of a equilateral H_4 ring is internally stable, but externally unstable. (As confirmed by the implementation of @xdvriend.)
 *  This test checks whether the RHF stability analysis confirms this.
 */
BOOST_AUTO_TEST_CASE(h4_sto3g_stability) {

    // Do our own RHF calculation
    const auto molecule = GQCP::Molecule::HRingFromDistance(4, 1.0);  // H3-triangle, 1 bohr apart

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "6-31G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in an AO basis
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian, 1.0e-05};

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain(1.0e-06, 1000);
    const auto qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    auto rhf_parameters = qc_structure.groundStateParameters();

    // We can now check the stability of the ground state parameters.
    // For this we need a restricted Hamiltonian in the orthonormal MO basis.
    const auto hamiltonian_restricted = sq_hamiltonian.transformed(rhf_parameters.expansion());

    // Calculate the stability matrices.
    const auto stability_matrices = rhf_parameters.calculateStabilityMatrices(hamiltonian_restricted);

    // This method should be internally stable.
    const auto internal_stability = stability_matrices.isInternallyStable();
    BOOST_CHECK(internal_stability == true);

    // This wavefunction should be externally unstable.
    const auto external_stability = stability_matrices.isExternallyStable();
    BOOST_CHECK(external_stability == false);

    // Both external stability subcases are checked individually as well.
    BOOST_CHECK(stability_matrices.isTripletStable() == false);
    BOOST_CHECK(stability_matrices.isComplexConjugateStable() == false);

    // Check that the stability properties can be printed
    stability_matrices.printStabilityDescription();
}
