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

#include "Basis/transform.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/Properties/expectation_values.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"

#include <random>


BOOST_AUTO_TEST_CASE(constrained_CO_test) {

    // Create a Molecule and the corresponding HamiltonianParameters
    auto CO = GQCP::Molecule::ReadXYZ("data/CO_mulliken.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {CO, "STO-3G"};
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, CO);  // in an AO basis

    size_t K = sq_hamiltonian.dimension();
    size_t N = CO.numberOfElectrons();

    GQCP::OneDM<double> D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);

    // Initialize the reference data from:
    // "Self-consistent methods constrained to a fixed number of particles in a given fragment and its relation to the electronegativity equalization method"
    // Authors: Andrés Cedillo, Dimitri Van Neck, Patrick Bultinck
    // DOI 10.1007/s00214-012-1227-6
    Eigen::Matrix<double, 21, 3> CO_data;
    //          lambda  charge  energy
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

    // Pick a set of AO's (GTOs of the carbon atom)
    std::vector<size_t> gto_list {0, 1, 2, 3, 4};

    // Iterate over multipliers
    size_t index = 0;
    for (double lambda = -1.0; lambda < 1.01; lambda += 0.1) {

        // Calculate the Mulliken operator
        auto mulliken_operator = spinor_basis.calculateMullikenOperator(gto_list);  // in AO basis

        // Constrain the original Hamiltonian
        auto constrained_sq_hamiltonian = sq_hamiltonian - lambda * mulliken_operator;  // in AO basis

        // Create a DIIS RHF SCF solver and solve the SCF equations
        auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(CO.numberOfElectrons(), constrained_sq_hamiltonian, spinor_basis.overlap().parameters());
        auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
        const GQCP::DiagonalRHFFockMatrixObjective<double> objective {constrained_sq_hamiltonian};
        const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment);
        const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

        // Transform only the mulliken operator to the basis in which the RHF energies are calculated
        mulliken_operator.transform(rhf_parameters.coefficientMatrix());

        // Retrieve the RHF "energy"
        double expectation_value = rhf_qc_structure.groundStateEnergy();

        // Retrieve the expectation value of the Mulliken operator (aka the population)
        double mulliken_population = mulliken_operator.calculateExpectationValue(D)(0);

        // Retrieve the total energy by adding the lambda times the expectation value of the constraining operator
        const double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(CO).value();
        double total_energy = expectation_value + lambda * mulliken_population + internuclear_repulsion_energy;

        // Mulliken charge on the carbon atom
        double C_charge = 6 - mulliken_population;

        // The lenient threshold is because of we are not certain of exact the geometry in the paper
        BOOST_CHECK(std::abs(total_energy - CO_data(index, 2)) < 1.0e-4);
        BOOST_CHECK(std::abs(C_charge - CO_data(index, 1)) < 1.0e-2);

        index++;
    }
}


BOOST_AUTO_TEST_CASE(constrained_CO_test_random_transformation) {
    // Repeat the same test but perform a random transformation
    // The Hartree-Fock procedure should be invariant under random transformations
    // This tests if our Mulliken operator remains correct if we transform before the procedure.

    // Create a Molecule and the corresponding HamiltonianParameters
    auto CO = GQCP::Molecule::ReadXYZ("data/CO_mulliken.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {CO, "STO-3G"};
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, CO);  // in an AO basis

    size_t K = sq_hamiltonian.dimension();
    size_t N = CO.numberOfElectrons();

    GQCP::TransformationMatrix<double> T = GQCP::TransformationMatrix<double>::Random(K, K);
    // set diagonal elements to 1
    for (int i = 0; i < K; i++) {
        T(i, i) = 1;
    }

    basisTransform(spinor_basis, sq_hamiltonian, T);

    GQCP::OneDM<double> D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);

    // Initialize the reference data from:
    // "Self-consistent methods constrained to a fixed number of particles in a given fragment and its relation to the electronegativity equalization method"
    // Authors: Andrés Cedillo, Dimitri Van Neck, Patrick Bultinck
    // DOI 10.1007/s00214-012-1227-6
    Eigen::Matrix<double, 21, 3> CO_data;
    //          lambda  charge  energy
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

    // Pick a set of AO's (GTOs of the carbon atom)
    std::vector<size_t> gto_list = {0, 1, 2, 3, 4};

    // Iterate over multipliers
    size_t index = 0;
    for (double lambda = -1; lambda < 1.01; lambda += 0.1) {

        // Calculate the Mulliken operator
        auto mulliken_operator = spinor_basis.calculateMullikenOperator(gto_list);

        // Contrain the original Hamiltonian
        auto sq_hamiltonian_constrained = sq_hamiltonian - lambda * mulliken_operator;

        // Create a DIIS RHF SCF solver and solve the SCF equations
        auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(CO.numberOfElectrons(), sq_hamiltonian_constrained, spinor_basis.overlap().parameters());
        auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
        const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian_constrained};
        const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment);
        const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

        // Transform only the Mulliken operator to the basis in which the RHF energies are calculated
        mulliken_operator.transform(rhf_parameters.coefficientMatrix());

        // Retrieve the RHF "energy"
        double expectation_value = rhf_qc_structure.groundStateEnergy();

        // Retrieve the expectation value of the Mulliken operator (aka the population)
        double mulliken_population = mulliken_operator.calculateExpectationValue(D)(0);

        // Retrieve the total energy by adding the lambda times the expectation value of the constraining operator
        double total_energy = expectation_value + lambda * mulliken_population + GQCP::Operator::NuclearRepulsion(CO).value();

        // Mulliken charge on the carbon atom
        double C_charge = 6 - mulliken_population;

        // The lenient threshold is because of we are not certain of exact the geometry in the paper
        BOOST_CHECK(std::abs(total_energy - CO_data(index, 2)) < 1.0e-4);
        BOOST_CHECK(std::abs(C_charge - CO_data(index, 1)) < 1.0e-2);

        index++;
    }
}
