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

#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/GHF/GHF.hpp"
#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"


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
    const auto S = g_spinor_basis.overlap();

    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


    // Create a solver and associated environment and let the QCMethod do its job.
    GQCP::SquareMatrix<double> C_initial_matrix {6};
    // clang-format off
    C_initial_matrix << -0.3585282,  0.0,        0.89935394,  0.0,         0.0,        1.57117404,
                        -0.3585282,  0.0,       -1.81035361,  0.0,         0.0,        0.00672366,
                        -0.3585282,  0.0,        0.91099966,  0.0,         0.0,        1.56445038,
                         0.0,       -0.3585282,  0.0,         0.89935394, -1.57117404, 0.0,
                         0.0,       -0.3585282,  0.0,        -1.81035361,  0.00672366, 0.0,
                         0.0,       -0.3585282,  0.0,         0.91099966,  1.56445038, 0.0;
    // clang-format on
    const GQCP::GTransformation<double> C_initial {C_initial_matrix};
    GQCP::GHFSCFEnvironment<double> environment {N, sq_hamiltonian, S, C_initial};

    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    const auto& ghf_parameters = qc_structure.groundStateParameters();


    // Provide reference values (from @xdvriend implementation) and check the results.
    const double ref_total_energy = -0.6311463202867755;
    GQCP::VectorX<double> ref_orbital_energies {6};
    ref_orbital_energies << -1.03323449, -0.89036198, 0.18717436, 0.76901002, 0.82078082, 0.93703294;

    const auto total_energy = qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-08);

    const auto orbital_energies = ghf_parameters.orbitalEnergies();
    BOOST_CHECK(orbital_energies.isApprox(ref_orbital_energies, 1.0e-06));


    // Check the reference value for S_x/y/z based on two different implementations.
    const double reference_s_x = 0.0;  // an UHF solution
    const double reference_s_y = 0.0;  // An UHF solution.
    const double reference_s_z = 0.5;  // An UHF solution.

    const GQCP::G1DM<GQCP::complex> P {ghf_parameters.calculateScalarBasis1DM().matrix().cast<GQCP::complex>()};  // The AO density matrix, converted from real to complex values.

    // Set up a complex spinor basis that can quantize the electronic spin operator.
    const GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> complex_g_spinor_basis {molecule, "STO-3G"};
    const auto S_op = complex_g_spinor_basis.quantize(GQCP::ElectronicSpinOperator());  // AO representation of the spin operator.

    const auto s_x = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::x);
    const auto s_y = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::y);
    const auto s_z1 = S_op.calculateExpectationValue(P)(GQCP::CartesianDirection::z);
    const auto s_z2 = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::z);

    BOOST_CHECK(std::abs(s_x - reference_s_x) < 1.0e-08);
    BOOST_CHECK(std::abs(s_y - reference_s_x) < 1.0e-08);
    BOOST_CHECK(std::abs(s_z1 - reference_s_z) < 1.0e-08);
    BOOST_CHECK(std::abs(s_z2 - reference_s_z) < 1.0e-08);
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
    const GQCP::GTransformation<double> C_initial {C_initial_matrix};
    GQCP::GHFSCFEnvironment<double> environment {N, sq_hamiltonian, S, C_initial};

    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 4000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    const auto& ghf_parameters = qc_structure.groundStateParameters();


    // Provide reference values (from @xdvriend implementation) and check the results.
    const double ref_total_energy = -0.6318365550450893;
    GQCP::VectorX<double> ref_orbital_energies {6};
    ref_orbital_energies << -0.96265264, -0.96265, 0.18610318, 0.82147981, 0.82147982, 0.8866645;

    const auto total_energy = qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-08);


    const auto orbital_energies = ghf_parameters.orbitalEnergies();
    BOOST_CHECK(orbital_energies.isApprox(ref_orbital_energies, 1.0e-06));


    // Check the reference value for S_x/y/z based on two different implementations.
    const double reference_s_x = 0.0;                     // a true GHF solution
    const double reference_s_y = 0.0;                     // A true GHF solution.
    const double reference_s_z = -4.903573113845816e-05;  // A true GHF solution.

    const GQCP::G1DM<GQCP::complex> P {ghf_parameters.calculateScalarBasis1DM().matrix().cast<GQCP::complex>()};  // The AO density matrix, converted from real to complex values.

    // Set up a complex spinor basis that can quantize the electronic spin operator.
    const GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> complex_g_spinor_basis {molecule, "STO-3G"};
    const auto S_op = complex_g_spinor_basis.quantize(GQCP::ElectronicSpinOperator());  // AO representation of the spin operator.

    const auto s_x = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::x);
    const auto s_y = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::y);
    const auto s_z1 = S_op.calculateExpectationValue(P)(GQCP::CartesianDirection::z);
    const auto s_z2 = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::z);

    BOOST_CHECK(std::abs(s_x - reference_s_x) < 1.0e-06);
    BOOST_CHECK(std::abs(s_y - reference_s_y) < 1.0e-06);
    BOOST_CHECK(std::abs(s_z1 - reference_s_z) < 1.0e-04);  // since the reference value is about 1.0e-05, this is the minimum threshold we can use
    BOOST_CHECK(std::abs(s_z2 - reference_s_z) < 1.0e-04);  // since the reference value is about 1.0e-05, this is the minimum threshold we can use
}


/**
 *  Check if the DIIS GHF SCF solver finds a correct solution.
 *
 *  The system of interest is a H3-triangle, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_test_DIIS) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.0);  // H3-triangle, 1 bohr apart
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


    // Create a solver and associated environment and let the QCMethod do its job.
    GQCP::SquareMatrix<double> C_initial_matrix {6};
    // clang-format off
    C_initial_matrix << -0.3585282,  0.0,        0.89935394,  0.0,         0.0,        1.57117404,
                        -0.3585282,  0.0,       -1.81035361,  0.0,         0.0,        0.00672366,
                        -0.3585282,  0.0,        0.91099966,  0.0,         0.0,        1.56445038,
                         0.0,       -0.3585282,  0.0,         0.89935394, -1.57117404, 0.0,
                         0.0,       -0.3585282,  0.0,        -1.81035361,  0.00672366, 0.0,
                         0.0,       -0.3585282,  0.0,         0.91099966,  1.56445038, 0.0;
    // clang-format on
    const GQCP::GTransformation<double> C_initial {C_initial_matrix};
    GQCP::GHFSCFEnvironment<double> environment {N, sq_hamiltonian, S, C_initial};

    auto solver = GQCP::GHFSCFSolver<double>::DIIS(6, 6, 1.0e-06, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    const auto ghf_parameters = qc_structure.groundStateParameters();


    // Provide reference values (from @xdvriend implementation) and check the results.
    const double ref_total_energy = -0.630521948908159;
    GQCP::VectorX<double> ref_orbital_energies {6};
    ref_orbital_energies << -1.03313925, -0.88946247, 0.18899685, 0.76709853, 0.81828059, 0.93860157;

    const auto total_energy = qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-08);

    const auto orbital_energies = ghf_parameters.orbitalEnergies();
    BOOST_CHECK(orbital_energies.isApprox(ref_orbital_energies, 1.0e-06));


    // Check the reference value for S_x/y/z based on two different implementations.
    const double reference_s_x = 0.0;                 // an UHF solution
    const double reference_s_y = 0.0;                 // An UHF solution.
    const double reference_s_z = 0.4999999999999999;  // An UHF solution.

    const GQCP::G1DM<GQCP::complex> P {ghf_parameters.calculateScalarBasis1DM().matrix().cast<GQCP::complex>()};  // The AO density matrix, converted from real to complex values.

    // Set up a complex spinor basis that can quantize the electronic spin operator.
    const GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> complex_g_spinor_basis {molecule, "STO-3G"};
    const auto S_op = complex_g_spinor_basis.quantize(GQCP::ElectronicSpinOperator());  // AO representation of the spin operator.

    const auto s_x = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::x);
    const auto s_y = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::y);
    const auto s_z1 = S_op.calculateExpectationValue(P)(GQCP::CartesianDirection::z);
    const auto s_z2 = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinOperator(), S)(GQCP::CartesianDirection::z);

    BOOST_CHECK(std::abs(s_x - reference_s_x) < 1.0e-08);
    BOOST_CHECK(std::abs(s_y - reference_s_x) < 1.0e-08);
    BOOST_CHECK(std::abs(s_z1 - reference_s_z) < 1.0e-08);
    BOOST_CHECK(std::abs(s_z2 - reference_s_z) < 1.0e-08);
}


/**
 *  Check if the specialized implementation of the GHF S^2 expectation value matches the general one.
 */
BOOST_AUTO_TEST_CASE(GHF_spin_expectation_values) {

    // Set up the molecular Hamiltonian in the AO basis.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.0);  // H3-triangle, 1 bohr apart.
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    const auto hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


    // Solve the GHF SCF equations, using a special initial guess from @xdvriend's implementation. This makes sure the off-diagonal spin-blocks are not zero blocks and helps the calculation to converge to a true GHF solution.
    GQCP::SquareMatrix<double> C_initial_matrix {6};
    // clang-format off
    C_initial_matrix << -0.3100721,  -0.15761163, -0.51612194, -0.38100148,  0.57090929, -0.37620802,
                        -0.00741269,  0.38801568, -0.25974834, -0.41043789, -0.67141074, -0.40332126,
                        -0.61961507,  0.18043708,  0.58367365,  0.17317687,  0.05464039, -0.45811451,
                         0.67031756,  0.28266352,  0.37079814, -0.23639173,  0.37758712, -0.3671939,
                         0.18059725, -0.8326703,   0.16282789, -0.03436191, -0.27832567, -0.41095738,
                         0.19477298,  0.13713633, -0.4018331,   0.77416187,  0.01572939, -0.42686445;
    // clang-format on
    const GQCP::GTransformation<double> C_initial {C_initial_matrix};
    GQCP::GHFSCFEnvironment<double> environment {N, hamiltonian, S, C_initial};

    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 4000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    const auto& ghf_parameters = qc_structure.groundStateParameters();


    // Provide reference values (from @xdvriend's implementation) and check the results.
    const double ref_total_energy = -0.6318365550450893;
    GQCP::VectorX<double> ref_orbital_energies {6};
    ref_orbital_energies << -0.96265264, -0.96265, 0.18610318, 0.82147981, 0.82147982, 0.8866645;

    const auto total_energy = qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-08);


    const auto orbital_energies = ghf_parameters.orbitalEnergies();
    BOOST_CHECK(orbital_energies.isApprox(ref_orbital_energies, 1.0e-06));


    // Calculate <S^2> through a specialized implementation.
    const auto s2_specialized = ghf_parameters.calculateExpectationValueOf(GQCP::ElectronicSpinSquaredOperator(), S);


    // Calculate <S^2> through a general expectation value, in MO basis. This requires the quantization of the S^2 operator, through some complex intermediates.
    const auto D_MO_real = ghf_parameters.calculateOrthonormalBasis1DM();
    const GQCP::G1DM<GQCP::complex> D_MO {D_MO_real.matrix().cast<GQCP::complex>()};

    auto d_MO_real = ghf_parameters.calculateOrthonormalBasis2DM();
    const GQCP::G2DM<GQCP::complex> d_MO {d_MO_real.tensor().cast<GQCP::complex>()};

    GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    GQCP::GTransformation<GQCP::complex> C {ghf_parameters.expansion().matrix().cast<GQCP::complex>()};
    spinor_basis.transform(C);
    const auto S2_MO = spinor_basis.quantize(GQCP::ElectronicSpinSquaredOperator());  // S^2 expressed in the GHF MO basis.

    const auto s2_general_MO = S2_MO.calculateExpectationValue(D_MO, d_MO);
    BOOST_CHECK(std::abs(s2_specialized - s2_general_MO) < 1.0e-08);


    // Calculate <S^2> through a general expectation value, in AO basis.
    // Note that the correct procedure to obtain S^2 in AO basis is to quantize S^2 in an orthonormal basis, and then back-transform it into AO basis.
    const auto D_AO = D_MO.transformed(C.inverse());
    const auto d_AO = d_MO.transformed(C.inverse());

    const auto S2_AO = S2_MO.transformed(C.inverse());

    const auto s2_general_AO = S2_AO.calculateExpectationValue(D_AO, d_AO);
    BOOST_CHECK(std::abs(s2_specialized - s2_general_AO) < 1.0e-08);
}

/**
 *  This test checks whether the lower lying complex GHF solution can indeed be found. 
 *  Note that this solution can also be found using real valued parameters.
 */
BOOST_AUTO_TEST_CASE(h3_sto3g_complex) {

    // Do our own GHF calculation.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.8897259886);  // H3-ring, 1 Angstrom apart.

    const GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> spinor_basis {molecule, "sto-3g"};
    const auto sq_hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In an AO basis.

    auto ghf_environment = GQCP::GHFSCFEnvironment<GQCP::complex>::WithComplexlyTransformedCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap());
    auto plain_ghf_scf_solver = GQCP::GHFSCFSolver<GQCP::complex>::Plain(1.0e-04, 1000);
    const auto qc_structure = GQCP::QCMethod::GHF<GQCP::complex>().optimize(plain_ghf_scf_solver, ghf_environment);
    auto ghf_parameters = qc_structure.groundStateParameters();
    auto ghf_ground_state_energy = qc_structure.groundStateEnergy();
    auto nuc_rep = GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();

    // Initialize a reference energy.
    const double reference_energy = -1.34044;

    // Both external stability subcases are checked individually as well.
    BOOST_CHECK(std::abs((ghf_ground_state_energy + nuc_rep) - reference_energy) < 1.0e-08);
}
