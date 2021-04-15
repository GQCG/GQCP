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

#define BOOST_TEST_MODULE "properties"

#include <boost/test/unit_test.hpp>

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/Transformations/transform.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationSolver.hpp"
#include "Operator/FirstQuantized/NuclearDipoleOperator.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Physical/units.hpp"
#include "Processing/Properties/RHFElectricalResponseSolver.hpp"
#include "Processing/Properties/properties.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"

/**
 *  Check the calculation of the CO dipole moment from a CCCBDB reference value.
 */
BOOST_AUTO_TEST_CASE(dipole_CO_STO_3G) {

    // Initialize the molecule and molecular Hamiltonian for CO.
    const GQCP::Nucleus C {6, 0.0, 0.0, 0.0};
    const GQCP::Nucleus O {8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145)};  // From CCCBDB, STO-3G geometry.
    const GQCP::Molecule molecule {{C, O}};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In the AO basis.

    const auto K = spin_orbital_basis.numberOfSpatialOrbitals();
    const size_t N = molecule.numberOfElectrons();

    // Solve the RHF SCF equations.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    const double total_energy = rhf_qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_REQUIRE(std::abs(total_energy - (-111.225)) < 1.0e-02);  // From CCCBDB, require a correct RHF solution to be found.


    // Calculate the RHF 1-DM and the dipole operator in RHF MO basis.
    const auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);
    auto dipole_op = spin_orbital_basis.quantize(GQCP::ElectronicDipoleOperator());
    dipole_op.transform(rhf_parameters.expansion());

    // Calculate the RHF total dipole moment in the MO basis and check with the reference value.
    GQCP::Vector<double, 3> total_dipole_moment = GQCP::NuclearDipoleOperator(molecule.nuclearFramework()).value() + dipole_op.calculateExpectationValue(D).asVector();
    BOOST_CHECK(std::abs(total_dipole_moment.norm() - (0.049)) < 1.0e-03);
}


/**
 *  Check the the RHF dipole moment for N2 is zero.
 */
BOOST_AUTO_TEST_CASE(dipole_N2_STO_3G) {

    // Initialize the molecule and the molecular Hamiltonian.
    const GQCP::Nucleus N_1 {7, 0.0, 0.0, 0.0};
    const GQCP::Nucleus N_2 {7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134)};  // From CCCBDB, STO-3G geometry.
    const GQCP::Molecule molecule {{N_1, N_2}};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In the AO basis.

    const auto K = spin_orbital_basis.numberOfSpatialOrbitals();
    const auto N = molecule.numberOfElectrons();

    // Solve the RHF SCF equations.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    const double total_energy = rhf_qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_REQUIRE(std::abs(total_energy - (-107.500654)) < 1.0e-05);  // From CCCBDB, require a correct RHF solution to be found.


    // Calculate the RHF 1-DM and the dipole operator in RHF MO basis.
    const auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);
    auto dipole_op = spin_orbital_basis.quantize(GQCP::ElectronicDipoleOperator());
    dipole_op.transform(rhf_parameters.expansion());

    // Calculate the RHF total dipole moment in the MO basis and check with the reference value.
    GQCP::Vector<double, 3> total_dipole_moment = GQCP::NuclearDipoleOperator(molecule.nuclearFramework()).value() + dipole_op.calculateExpectationValue(D).asVector();
    BOOST_CHECK(std::abs(total_dipole_moment.norm() - (0.0)) < 1.0e-08);
}


/**
 *  Check the calculation of the zz-component of the polarizability for H2 with a reference value from Psi4-numpy.
 * 
 *  Note that the reference value is generated from Psi4-numpy, with a fix for the Fockian matrix.
 */
BOOST_AUTO_TEST_CASE(h2_polarizability_RHF) {

    // Initialize the reference value.
    const double ref_alpha_zz = 1.08428;


    // Initialize the molecule and the Hamiltonian in the AO basis.
    const GQCP::Nucleus H1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus H2 {1, 0.0, 0.0, 0.5};
    const GQCP::Molecule molecule {{H1, H2}, 0};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In the AO basis.


    // Do the RHF calculation to get the canonical RHF orbitals.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective(sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();


    // Transform the orbitals to the RHF basis and prepare the dipole integrals in the RHF basis.
    GQCP::transform(rhf_parameters.expansion(), spin_orbital_basis, sq_hamiltonian);
    const auto dipole_op = spin_orbital_basis.quantize(GQCP::ElectronicDipoleOperator());


    // Find the RHF wave function response.
    GQCP::RHFElectricalResponseSolver cphf_solver {molecule.numberOfElectrons() / 2};
    const auto x = cphf_solver.calculateWaveFunctionResponse(sq_hamiltonian, dipole_op);


    // Calculate the RHF polarizability and check with the reference value.
    const auto F_p = cphf_solver.calculateParameterResponseForce(dipole_op);
    const auto alpha = GQCP::calculateElectricPolarizability(F_p, x);
    const auto alpha_zz = alpha(2, 2);

    BOOST_CHECK(std::abs(alpha_zz - ref_alpha_zz) < 1.0e-05);
}


/**
 *  Check the calculation of the ipsocentric current density and the intermediates for its calculation. The test system is H2, 1 au apart in an STO-3G basis set.
 *
 *  The reference implementation is a GAMESS-SYSMO combination. Calculations were performed by Remco Havenith.
 */
BOOST_AUTO_TEST_CASE(ipsocentric_H2) {

    using namespace GQCP::literals;

    // Set up the molecular Hamiltonian in AO basis.
    const GQCP::Molecule molecule {{GQCP::Nucleus(1, 0.0, 0.0, 0.0), GQCP::Nucleus(1, 0.0, 0.0, 1.0)}};
    const auto N = molecule.numberOfElectrons();
    const auto N_P = molecule.numberOfElectronPairs();

    const std::string basis_set {"STO-3G"};
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, basis_set};

    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


    // Read in the GAMESS-UK RHF wave function model parameters. Even though GQCP finds the same orbitals and orbital energies, the phase factors (and orbital coefficients for degenerated orbitals) are unlikely to be reproducible.
    GQCP::MatrixX<double> C_matrix {2, 2};
    // clang-format off
    C_matrix << 0.5275464665,  1.5678230259,
                0.5275464665, -1.5678230259;
    // clang-format on
    GQCP::RTransformation<double> C {C_matrix};

    GQCP::VectorX<double> orbital_energies {2};
    orbital_energies << -0.6757801904, 0.9418115528;

    const GQCP::QCModel::RHF<double> rhf_parameters {N_P, orbital_energies, C};


    // Since we're going to work with complex operators, we have to let a complex spin-orbital basis do the quantization.
    GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::GTOShell> complex_spin_orbital_basis {molecule, basis_set};
    GQCP::RTransformation<GQCP::complex> C_complex {C_matrix.cast<GQCP::complex>()};
    complex_spin_orbital_basis.transform(C_complex);

    spin_orbital_basis.transform(C);
    hamiltonian.transform(C);


    // Calculate the orbital Hessian.
    const auto orbital_space = rhf_parameters.orbitalSpace();
    auto A = rhf_parameters.calculateOrbitalHessianForImaginaryResponse(hamiltonian, orbital_space);

    GQCP::MatrixX<double> A_ref {1, 1};
    // clang-format off
    A_ref << 2.17196;
    // clang-format on

    BOOST_CHECK(A_ref.isApprox((A.asMatrix() * (0.5_ii)).real(), 1.0e-04));


    // Solve the CPHF equations for the angular momentum operator.
    const auto L = complex_spin_orbital_basis.quantize(GQCP::AngularMomentumOperator());
    const auto F_B = rhf_parameters.calculateMagneticFieldResponseForce(L);

    GQCP::MatrixX<double> F_B_ref {3, 1};
    // clang-format off
    F_B_ref << 0.00000,
               0.00000,
               0.00000;
    // clang-format on

    BOOST_CHECK(F_B_ref.transpose().isApprox((F_B * (0.5_ii)).real(), 1.0e-04));


    auto environment_B = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_B);
    auto solver_B = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_B.perform(environment_B);
    const auto x = environment_B.x;

    GQCP::MatrixX<double> x_ref {3, 1};
    // clang-format off
    x_ref << 0.00000,
             0.00000,
             0.00000;
    // clang-format on

    BOOST_CHECK(x_ref.transpose().isApprox(-x.real(), 1.0e-04));


    // Solve the CPHF equations for the linear momentum operator.
    const auto p = complex_spin_orbital_basis.quantize(GQCP::LinearMomentumOperator());
    const auto F_G = rhf_parameters.calculateGaugeOriginTranslationResponseForce(p);

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
    GQCP::MatrixX<GQCP::complex> F_G_reduced {1, 3};
    F_G_reduced.col(0) = 2 * F_G.col(3);  // x <--> yz
    F_G_reduced.col(1) = 2 * F_G.col(4);  // y <--> zx
    F_G_reduced.col(2) = 2 * F_G.col(0);  // z <--> xy

    GQCP::MatrixX<double> F_G_ref {3, 1};
    // clang-format off
    F_G_ref <<  0.00000,
                0.00000,
               -1.09788;
    // clang-format on

    BOOST_CHECK(F_G_ref.transpose().isApprox((F_G_reduced * (0.5_ii)).real(), 1.0e-04));


    auto environment_G = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_G);
    auto solver_G = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_G.perform(environment_G);

    const auto y = environment_G.x;

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
    GQCP::MatrixX<GQCP::complex> y_reduced {1, 3};
    y_reduced.col(0) = 2 * y.col(3);  // x <--> yz
    y_reduced.col(1) = 2 * y.col(4);  // y <--> zx
    y_reduced.col(2) = 2 * y.col(0);  // z <--> xy

    GQCP::MatrixX<double> y_ref {3, 1};
    // clang-format off
    y_ref <<  0.00000,
              0.00000,
             -0.50548;
    // clang-format on

    BOOST_CHECK(y_ref.transpose().isApprox(-y_reduced.real(), 1.0e-04));


    // Calculate the ipsocentric magnetic inducibility on a grid and check the results.
    const auto j_op = complex_spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());

    const GQCP::Vector<double, 3> origin {-2.0, -2.0, -2.0};
    const std::array<size_t, 3> steps {6, 6, 6};
    const std::array<double, 3> step_sizes {0.8, 0.8, 0.8};

    const GQCP::CubicGrid grid {origin, steps, step_sizes};
    const auto grid_points = grid.points();

    const auto J_field = GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(grid, orbital_space, x, y, j_op);
    const auto& J_field_values = J_field.values();

    const auto J_field_ref = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/h2_sysmo_currents.rgrid");
    const auto& J_field_ref_values = J_field_ref.values();
    for (size_t i = 0; i < grid.numberOfPoints(); i++) {
        std::cout << grid_points[i] << std::endl;
        std::cout << "----------" << std::endl;

        std::cout << J_field_values[i].col(GQCP::CartesianDirection::y).real() << std::endl << std::endl;
        
        std::cout << J_field_ref_values[i] << std::endl << std::endl;

        // BOOST_CHECK(J_field_values[i].col(GQCP::CartesianDirection::y).real().isApprox(J_field_ref_values[i], 1.0e-08));  // The reference values contain the magnetic field response in the y-direction?
    }
}


// /**
//  *  Check the calculation of the ipsocentric current density and the intermediates for its calculation. The test system is H2O in an STO-3G basis set.
//  *
//  *  The reference implementation is a GAMESS-SYSMO combination. Calculations were performed by Remco Havenith.
//  */
// // BOOST_AUTO_TEST_CASE(ipsocentric_H2O) {

// //     using namespace GQCP::literals;

// //     // Set up the molecular Hamiltonian in AO basis.
// //     const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
// //     const auto N = molecule.numberOfElectrons();
// //     const auto N_P = molecule.numberOfElectronPairs();

// //     const std::string basis_set {"STO-3G"};
// //     GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, basis_set};

// //     auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


// //     // Solve the RHF SCF equations.
// //     const auto S = spin_orbital_basis.quantize(GQCP::OverlapOperator());
// //     auto scf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, hamiltonian, S);
// //     auto scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
// //     const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};

// //     const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, scf_solver, scf_environment).groundStateParameters();


// //     // Adjust the RHF model parameters to match the solution found by GAMESS. Even though GQCP finds the same orbitals, up to a phase factor, we have to adjust them to be able to check further calculations.
// //     const auto C = rhf_parameters.expansion();

// //     auto C_adjusted_matrix = C.matrix();
// //     C_adjusted_matrix.col(1) *= -1;
// //     C_adjusted_matrix.col(3) *= -1;
// //     C_adjusted_matrix.col(4) *= -1;
// //     C_adjusted_matrix.col(5) *= -1;
// //     C_adjusted_matrix.col(6) *= -1;
// //     GQCP::RTransformation<double> C_adjusted {C_adjusted_matrix};

// //     GQCP::MatrixX<double> C_adjusted_ref {7, 7};
// //     // clang-format off
// //     C_adjusted_ref << -0.9944,  0.2392,  0.0000,  0.0937, 0.0000, -0.1116,  0.0000,
// //                       -0.0241, -0.8857,  0.0000, -0.4796, 0.0000,  0.6696, -0.0000,
// //                       -0.0000, -0.0000, -0.6073, -0.0000, 0.0000,  0.0000,  0.9192,
// //                       -0.0032, -0.0859, -0.0000,  0.7474, 0.0000,  0.7385, -0.0000,
// //                        0.0000,  0.0000,  0.0000,  0.0000, 1.0000,  0.0000,  0.0000,
// //                        0.0046, -0.1440, -0.4530,  0.3295, 0.0000, -0.7098, -0.7325,
// //                        0.0046, -0.1440,  0.4530,  0.3295, 0.0000, -0.7098,  0.7325;
// //     // clang-format on
// //     BOOST_REQUIRE(C_adjusted_ref.isApprox(C_adjusted.matrix(), 1.0e-03));

// //     const GQCP::QCModel::RHF<double> rhf_parameters_adjusted {N_P, rhf_parameters.orbitalEnergies(), C_adjusted};


// //     // Since we're going to work with complex operators, we have to let a complex spin-orbital basis do the quantization.
// //     GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::GTOShell> complex_spin_orbital_basis {molecule, basis_set};
// //     GQCP::RTransformation<GQCP::complex> C_adjusted_complex {C_adjusted.matrix().cast<GQCP::complex>()};
// //     complex_spin_orbital_basis.transform(C_adjusted_complex);

// //     spin_orbital_basis.transform(C_adjusted);
// //     hamiltonian.transform(C_adjusted);


// //     // Calculate the orbital Hessian.
// //     const auto orbital_space = rhf_parameters.orbitalSpace();
// //     auto A = rhf_parameters_adjusted.calculateOrbitalHessianForImaginaryResponse(hamiltonian, orbital_space);

// //     GQCP::MatrixX<double> A_ref {10, 10};
// //     // clang-format off
// //     A_ref << 39.96940,  0.00000,  0.03448,  0.00000,  0.00000, -0.00925,  0.03956,  0.00000,  0.00000,  0.00000,
// //               0.00000, 40.06171,  0.00000,  0.05605, -0.02686,  0.00000,  0.00000,  0.02998,  0.00000,  0.00000,
// //               0.03448,  0.00000,  2.37666,  0.00000,  0.00000, -0.11291, -0.02343,  0.00000,  0.00000,  0.00000,
// //               0.00000,  0.05605,  0.00000,  2.51331, -0.03566,  0.00000,  0.00000, -0.00228,  0.00000,  0.00000,
// //               0.00000, -0.02686,  0.00000, -0.03566,  1.10153,  0.00000,  0.00000, -0.07674,  0.00000,  0.00000,
// //              -0.00925,  0.00000, -0.11291,  0.00000,  0.00000,  1.42443,  0.01478,  0.00000,  0.00000,  0.00000,
// //               0.03956,  0.00000, -0.02343,  0.00000,  0.00000,  0.01478,  0.93180,  0.00000,  0.00000,  0.00000,
// //               0.00000,  0.02998,  0.00000, -0.00228, -0.07674,  0.00000,  0.00000,  1.03049,  0.00000,  0.00000,
// //               0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.64372,  0.00000,
// //               0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.78206;
// //     // clang-format on

// //     BOOST_CHECK(A_ref.isApprox((A.asMatrix() * (0.5_ii)).real(), 1.0e-04));


// //     // Solve the CPHF equations for the angular momentum operator.
// //     const auto L = complex_spin_orbital_basis.quantize(GQCP::AngularMomentumOperator());
// //     const auto F_B = rhf_parameters_adjusted.calculateMagneticFieldResponseForce(L);

// //     GQCP::MatrixX<double> F_B_ref {3, 10};
// //     // clang-format off
// //     F_B_ref << 0.00000, 0.00000, 0.00000,  0.00000, 0.00000,  0.00000, 0.00000, 0.00000,  0.42752, -0.00000,
// //                0.00000, 0.00000, 0.00000,  0.00000, 0.00000,  0.00000, 0.00000, 0.00000, -0.00000, -0.52599,
// //                0.00000, 0.12303, 0.00000, -0.06987, 0.26617, -0.00000, 0.00000, 0.47815,  0.00000,  0.00000;
// //     // clang-format on

// //     BOOST_CHECK(F_B_ref.transpose().isApprox((F_B * (0.5_ii)).real(), 1.0e-04));


// //     auto environment_B = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_B);
// //     auto solver_B = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
// //     solver_B.perform(environment_B);
// //     const auto x = environment_B.x;

// //     GQCP::MatrixX<double> x_ref {3, 10};
// //     // clang-format off
// //     x_ref << 0.00000, 0.00000, 0.00000,  0.00000, 0.00000,  0.00000, 0.00000, 0.00000,  0.66414, -0.00000,
// //              0.00000, 0.00000, 0.00000,  0.00000, 0.00000,  0.00000, 0.00000, 0.00000, -0.00000, -0.67257,
// //              0.00000, 0.00293, 0.00000, -0.02353, 0.27468, -0.00000, 0.00000, 0.48432,  0.00000,  0.00000;
// //     // clang-format on

// //     BOOST_CHECK(x_ref.transpose().isApprox(-x.real(), 1.0e-04));


// //     // Solve the CPHF equations for the linear momentum operator.
// //     const auto p = complex_spin_orbital_basis.quantize(GQCP::LinearMomentumOperator());
// //     const auto F_G = rhf_parameters_adjusted.calculateGaugeOriginTranslationResponseForce(p);

// //     // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
// //     GQCP::MatrixX<GQCP::complex> F_G_reduced {10, 3};
// //     F_G_reduced.col(0) = 2 * F_G.col(3);  // x <--> yz
// //     F_G_reduced.col(1) = 2 * F_G.col(4);  // y <--> zx
// //     F_G_reduced.col(2) = 2 * F_G.col(0);  // z <--> xy

// //     GQCP::MatrixX<double> F_G_ref {3, 10};
// //     // clang-format off
// //     F_G_ref << -0.00000, -1.73412, 0.00000, 0.12306, 0.84094, 0.00000, -0.00000, -0.67978,  0.00000, 0.00000,
// //                -1.39805,  0.00000, 0.06859, 0.00000, 0.00000, 0.61259, -0.64447,  0.00000,  0.00000, 0.00000,
// //                 0.00000,  0.00000, 0.00000, 0.00000, 0.00000, 0.00000,  0.00000,  0.00000, -0.18459, 0.00000;
// //     // clang-format on

// //     BOOST_CHECK(F_G_ref.transpose().isApprox((F_G_reduced * (0.5_ii)).real(), 1.0e-04));


// //     auto environment_G = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_G);
// //     auto solver_G = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
// //     solver_G.perform(environment_G);

// //     const auto y = environment_G.x;
// //     // std::cout << y << std::endl;

// //     // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
// //     GQCP::MatrixX<GQCP::complex> y_reduced {10, 3};
// //     y_reduced.col(0) = 2 * y.col(3);  // x <--> yz
// //     y_reduced.col(1) = 2 * y.col(4);  // y <--> zx
// //     y_reduced.col(2) = 2 * y.col(0);  // z <--> xy

// //     GQCP::MatrixX<double> y_ref {3, 10};
// //     // clang-format off
// //     y_ref << -0.00000, -0.04243, 0.00000, 0.05961, 0.72221, 0.00000, -0.00000, -0.60452,  0.00000, 0.00000,
// //              -0.03422,  0.00000, 0.04342, 0.00000, 0.00000, 0.44050, -0.69608,  0.00000,  0.00000, 0.00000,
// //               0.00000,  0.00000, 0.00000, 0.00000, 0.00000, 0.00000,  0.00000,  0.00000, -0.28675, 0.00000;
// //     // clang-format on

// //     BOOST_CHECK(y_ref.transpose().isApprox(-y_reduced.real(), 1.0e-04));


// //     const auto j_op = complex_spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());

// //     std::cout << "j_op calculated" << std::endl;

// //     // const GQCP::Vector<double, 3> origin {-2.0, -2.0, 0.0};
// //     // const std::array<size_t, 3> steps {24, 24, 1};
// //     // const std::array<double, 3> step_sizes {4.0 / 24, 4.0 / 24, 0.0};

// //     // const GQCP::CubicGrid grid {origin, steps, step_sizes};

// //     // const auto J_field = GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(grid, orbital_space, x, y, j_op);

// //     for (size_t i = 0; i < 25; i++) {
// //         GQCP::Vector<double, 3> r {-2.0 + i * 4.0 / 24, -2.0, 0.0};
// //         std::cout << r << std::endl;
// //         std::cout << "----------" << std::endl;
// //         std::cout << GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(GQCP::Vector<double, 3>(-2.0 + i * 4.0 / 24, -2.0, 0.0), orbital_space, x, y, j_op).col(2).real() << std::endl
// //                   << std::endl;
// //     }


// //     // const auto& values = J_field.values();
// //     // for (const auto& value : values) {
// //     //     std::cout << value.real().col(2) << std::endl
// //     //               << std::endl;
// //     // }
// // }
