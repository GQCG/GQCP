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

#define BOOST_TEST_MODULE "constrained_RHF"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/transform.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"

#include <random>


/**
 *  Check constrained RHF Mulliken populations, based on the reference article: "Self-consistent methods constrained to a fixed number of particles in a given fragment and its relation to the electronegativity equalization method", by Andrés Cedillo, Dimitri Van Neck, Patrick Bultinck. (DOI: 10.1007/s00214-012-1227-6).
 * 
 *  We still use a lenient threshold because we are not certain of exact the geometry in the paper.
 */
BOOST_AUTO_TEST_CASE(constrained_CO_test) {

    // The reference data from the article (see below). The first column represents the multiplier, the second column the Mulliken charge on 'C', the third column the constrained RHF energy.
    GQCP::Matrix<double, 21, 3> CO_data;
    // clang-format off
    CO_data <<  -1.0,  1.73, -110.530475,
                -0.9,  1.62, -110.634766,
                -0.8,  1.50, -110.737606,
                -0.7,  1.37, -110.836596,
                -0.6,  1.23, -110.929219,
                -0.5,  1.08, -111.012983,
                -0.4,  0.91, -111.085538,
                -0.3,  0.75, -111.144760,
                -0.2,  0.57, -111.188802,
                -0.1,  0.39, -111.216114,
                 0.0,  0.20, -111.225446,
                 0.1,  0.01, -111.215842,
                 0.2, -0.19, -111.186619,
                 0.3, -0.38, -111.137356,
                 0.4, -0.58, -111.067872,
                 0.5, -0.78, -110.978210,
                 0.6, -0.98, -110.868621,
                 0.7, -1.18, -110.739544,
                 0.8, -1.38, -110.591599,
                 0.9, -1.57, -110.425574,
                 1.0, -1.77, -110.242423;
    // clang-format on


    // Prepare the molecular Hamiltonian in the AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/CO_mulliken.xyz");
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    const auto S = spin_orbital_basis.overlap();

    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.
    const auto K = hamiltonian.numberOfOrbitals();


    // Prepare the Mulliken partitioning on 'C'.
    using Shell = GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell>::Shell;
    const auto mulliken_partitioning = spin_orbital_basis.mullikenPartitioning(
        [](const Shell& shell) {
            return shell.nucleus().element() == "C";
        });


    // Iterate over the Lagrange multipliers and check if we can reproduce all values.
    size_t index = 0;
    for (double lambda = -1.0; lambda < 1.01; lambda += 0.1) {

        // Do the constrained RHF calculation by adding (-lambda * Mulliken) to the molecular Hamiltonian.
        auto mulliken_op = S.partitioned(mulliken_partitioning);
        const auto constrained_hamiltonian = hamiltonian - lambda * mulliken_op;  // In the AO basis.

        // Create a DIIS RHF SCF solver and solve the SCF equations.
        auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), constrained_hamiltonian, S);
        auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
        const GQCP::DiagonalRHFFockMatrixObjective<double> objective {constrained_hamiltonian};

        const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment);
        const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

        // Calculate the Mulliken population through the expectation value of the Mulliken operator. We can perform this calculation in the AO basis.
        const auto D_RHF = rhf_parameters.calculateScalarBasis1DM();
        const double mulliken_population = mulliken_op.calculateExpectationValue(D_RHF);  // A scalar-StorageArray can be implicitly converted into the underlying scalar.


        // Calculate the total internuclear energy as a contribution from three sources, and check with the reference value.
        const double internuclear_repulsion_energy = GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
        double total_energy = rhf_qc_structure.groundStateEnergy() + lambda * mulliken_population + internuclear_repulsion_energy;

        BOOST_CHECK(std::abs(total_energy - CO_data(index, 2)) < 1.0e-4);


        // Calculate the Mulliken charge on 'C' and check with the reference value.
        const double C_charge = 6 - mulliken_population;
        BOOST_CHECK(std::abs(C_charge - CO_data(index, 1)) < 1.0e-2);

        index++;
    }
}


/**
 *  Repeat the previous test, but use a random non-orthogonal basis. It'll test the correct implementation of the Mulliken operator (aka the Mulliken-partitioned number operator).
 */
BOOST_AUTO_TEST_CASE(constrained_CO_test_random_AO_basis) {

    // The reference data from the article (see below). The first column represents the multiplier, the second column the Mulliken charge on 'C', the third column the constrained RHF energy.
    GQCP::Matrix<double, 21, 3> CO_data;
    // clang-format off
    CO_data <<  -1.0,  1.73, -110.530475,
                -0.9,  1.62, -110.634766,
                -0.8,  1.50, -110.737606,
                -0.7,  1.37, -110.836596,
                -0.6,  1.23, -110.929219,
                -0.5,  1.08, -111.012983,
                -0.4,  0.91, -111.085538,
                -0.3,  0.75, -111.144760,
                -0.2,  0.57, -111.188802,
                -0.1,  0.39, -111.216114,
                 0.0,  0.20, -111.225446,
                 0.1,  0.01, -111.215842,
                 0.2, -0.19, -111.186619,
                 0.3, -0.38, -111.137356,
                 0.4, -0.58, -111.067872,
                 0.5, -0.78, -110.978210,
                 0.6, -0.98, -110.868621,
                 0.7, -1.18, -110.739544,
                 0.8, -1.38, -110.591599,
                 0.9, -1.57, -110.425574,
                 1.0, -1.77, -110.242423;
    // clang-format on


    // Prepare the molecular Hamiltonian in a random AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/CO_mulliken.xyz");
    const auto N = molecule.numberOfElectrons();
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    const auto K = spin_orbital_basis.numberOfSpatialOrbitals();

    // Generate a random transformation, but keep the norm of the orbitals intact.
    GQCP::SquareMatrix<double> T_matrix = GQCP::SquareMatrix<double>::Random(K);
    for (size_t i = 0; i < K; i++) {
        T_matrix(i, i) = 1.0;
    }
    const GQCP::RTransformation<double> T {T_matrix};

    spin_orbital_basis.transform(T);
    const auto S = spin_orbital_basis.overlap();

    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.


    // Prepare the Mulliken partitioning on 'C'.
    using Shell = GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell>::Shell;
    const auto mulliken_partitioning = spin_orbital_basis.mullikenPartitioning(
        [](const Shell& shell) {
            return shell.nucleus().element() == "C";
        });


    // Iterate over the Lagrange multipliers and check if we can reproduce all values.
    size_t index = 0;
    for (double lambda = -1.0; lambda < 1.01; lambda += 0.1) {

        // Do the constrained RHF calculation by adding (-lambda * Mulliken) to the molecular Hamiltonian.
        auto mulliken_op = S.partitioned(mulliken_partitioning);
        const auto constrained_hamiltonian = hamiltonian - lambda * mulliken_op;  // In the AO basis.

        // Create a DIIS RHF SCF solver and solve the SCF equations.
        auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), constrained_hamiltonian, S);
        auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
        const GQCP::DiagonalRHFFockMatrixObjective<double> objective {constrained_hamiltonian};

        const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment);
        const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

        // Calculate the Mulliken population through the expectation value of the Mulliken operator. We can perform this calculation in the AO basis.
        const auto D_RHF = rhf_parameters.calculateScalarBasis1DM();
        const double mulliken_population = mulliken_op.calculateExpectationValue(D_RHF);  // A scalar-StorageArray can be implicitly converted into the underlying scalar.


        // Calculate the total internuclear energy as a contribution from three sources, and check with the reference value.
        const double internuclear_repulsion_energy = GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
        double total_energy = rhf_qc_structure.groundStateEnergy() + lambda * mulliken_population + internuclear_repulsion_energy;

        BOOST_CHECK(std::abs(total_energy - CO_data(index, 2)) < 1.0e-4);


        // Calculate the Mulliken charge on 'C' and check with the reference value.
        const double C_charge = 6 - mulliken_population;
        BOOST_CHECK(std::abs(C_charge - CO_data(index, 1)) < 1.0e-2);

        index++;
    }
}
