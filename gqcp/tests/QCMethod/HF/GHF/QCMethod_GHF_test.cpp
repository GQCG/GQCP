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

#define BOOST_TEST_MODULE "QCMethod_GHF_test"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/GHF/GHF.hpp"
#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"
#include "Utilities/linalg.hpp"


/**
 *  The following tests check if the default implementations for the GHF SCF solvers give the correct result.
 */

/**
 *  Check if the plain GHF SCF solver finds a correct solution.
 * 
 *  The system of interest is a H3-triangle, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_test_1) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.0);  // H3-triangle, 1 bohr apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap().parameters();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(g_spinor_basis, molecule);


    // Create a solver and associated environment and let the QCMethod do its job.
    GQCP::TransformationMatrix<double> C_initial {6};
    // clang-format off
    C_initial << -0.3585282,  0.0,        0.89935394,  0.0,         0.0,        1.57117404,
                 -0.3585282,  0.0,       -1.81035361,  0.0,         0.0,        0.00672366,
                 -0.3585282,  0.0,        0.91099966,  0.0,         0.0,        1.56445038,
                  0.0,       -0.3585282,  0.0,         0.89935394, -1.57117404, 0.0,
                  0.0,       -0.3585282,  0.0,        -1.81035361,  0.00672366, 0.0,
                  0.0,       -0.3585282,  0.0,         0.91099966,  1.56445038, 0.0;
    // clang-format on
    GQCP::GHFSCFEnvironment<double> environment {N, sq_hamiltonian, S, C_initial};

    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    const auto ghf_parameters = qc_structure.groundStateParameters();


    // Provide reference values (from @xdvriend implementation) and check the results.
    const double ref_total_energy = -0.6311463202867755;
    GQCP::VectorX<double> ref_orbital_energies {6};
    ref_orbital_energies << -1.03323449, -0.89036198, 0.18717436, 0.76901002, 0.82078082, 0.93703294;

    const auto total_energy = qc_structure.groundStateEnergy() + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-08);

    const auto orbital_energies = ghf_parameters.orbitalEnergies();
    BOOST_CHECK(orbital_energies.isApprox(ref_orbital_energies, 1.0e-06));


    // Check the reference value for S_z based on two different implementations.
    // const double reference_s_z = 0.5;  // an UHF solution

    // const auto P = ghf_parameters.calculateScalarBasis1RDM();                     // AO density matrix
    // const auto S_op = g_spinor_basis.quantize(GQCP::Operator::ElectronicSpin());  // AO representation of the spin operator

    // const auto s_z1 = S_op.calculateExpectationValue(P)(GQCP::CartesianDirection::z);
    // const auto s_z2 = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::z);

    // BOOST_CHECK(std::abs(s_z1 - reference_s_z) < 1.0e-08);
    // BOOST_CHECK(std::abs(s_z2 - reference_s_z) < 1.0e-08);
}


/**
 *  Check if the plain GHF SCF solver finds a correct solution.
 * 
 *  The system of interest is a H3-triangle, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_test_2) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.0);  // H3-triangle, 1 bohr apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap().parameters();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(g_spinor_basis, molecule);


    // Create a solver and associated environment and let the QCMethod do its job.
    GQCP::TransformationMatrix<double> C_initial {6};
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
    const auto ghf_parameters = qc_structure.groundStateParameters();


    // Provide reference values (from @xdvriend implementation) and check the results.
    const double ref_total_energy = -0.6318365550450893;
    GQCP::VectorX<double> ref_orbital_energies {6};
    ref_orbital_energies << -0.96265264, -0.96265, 0.18610318, 0.82147981, 0.82147982, 0.8866645;

    const auto total_energy = qc_structure.groundStateEnergy() + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-08);


    const auto orbital_energies = ghf_parameters.orbitalEnergies();
    BOOST_CHECK(orbital_energies.isApprox(ref_orbital_energies, 1.0e-06));


    // Check the reference value for S_z based on two different implementations.
    // const double reference_s_z = -4.903573113845816e-05;  // a true GHF solution

    // const auto P = ghf_parameters.calculateScalarBasis1RDM();                     // AO density matrix
    // const auto S_op = g_spinor_basis.quantize(GQCP::Operator::ElectronicSpin());  // AO representation of the spin operator

    // const auto s_z1 = S_op.calculateExpectationValue(P)(GQCP::CartesianDirection::z);
    // const auto s_z2 = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::z);

    // BOOST_CHECK(std::abs(s_z1 - reference_s_z) < 1.0e-04);  // since the reference value is about 1.0e-05, this is the minimum threshold we can use
    // BOOST_CHECK(std::abs(s_z2 - reference_s_z) < 1.0e-04);  // since the reference value is about 1.0e-05, this is the minimum threshold we can use
}


/**
 *  Check if the plain GHF SCF solver finds a correct solution.
 * 
 *  The system of interest is a H3-triangle, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_test_DIIS) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.0);  // H3-triangle, 1 bohr apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap().parameters();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(g_spinor_basis, molecule);


    // Create a solver and associated environment and let the QCMethod do its job.
    GQCP::TransformationMatrix<double> C_initial {6};
    // clang-format off
    C_initial << -0.3585282,  0.0,        0.89935394,  0.0,         0.0,        1.57117404,
                 -0.3585282,  0.0,       -1.81035361,  0.0,         0.0,        0.00672366,
                 -0.3585282,  0.0,        0.91099966,  0.0,         0.0,        1.56445038,
                  0.0,       -0.3585282,  0.0,         0.89935394, -1.57117404, 0.0,
                  0.0,       -0.3585282,  0.0,        -1.81035361,  0.00672366, 0.0,
                  0.0,       -0.3585282,  0.0,         0.91099966,  1.56445038, 0.0;
    // clang-format on
    GQCP::GHFSCFEnvironment<double> environment {N, sq_hamiltonian, S, C_initial};

    auto solver = GQCP::GHFSCFSolver<double>::DIIS(6, 6, 1.0e-06, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    const auto ghf_parameters = qc_structure.groundStateParameters();


    // Provide reference values (from @xdvriend implementation) and check the results.
    const double ref_total_energy = -0.630521948908159;
    GQCP::VectorX<double> ref_orbital_energies {6};
    ref_orbital_energies << -1.03313925, -0.88946247, 0.18899685, 0.76709853, 0.81828059, 0.93860157;

    const auto total_energy = qc_structure.groundStateEnergy() + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-08);

    const auto orbital_energies = ghf_parameters.orbitalEnergies();
    BOOST_CHECK(orbital_energies.isApprox(ref_orbital_energies, 1.0e-06));


    // Check the reference value for S_z based on two different implementations.
    // const double reference_s_z = 0.4999999999999999;  // an UHF solution

    // const auto P = ghf_parameters.calculateScalarBasis1RDM();                     // AO density matrix
    // const auto S_op = g_spinor_basis.quantize(GQCP::Operator::ElectronicSpin());  // AO representation of the spin operator

    // const auto s_z1 = S_op.calculateExpectationValue(P)(GQCP::CartesianDirection::z);
    // const auto s_z2 = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::z);

    // BOOST_CHECK(std::abs(s_z1 - reference_s_z) < 1.0e-08);
    // BOOST_CHECK(std::abs(s_z2 - reference_s_z) < 1.0e-08);
}
