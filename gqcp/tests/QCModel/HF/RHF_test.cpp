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

#define BOOST_TEST_MODULE "QCModel::RHF"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCModel/HF/RHF.hpp"


/**
 *  Check for valid and invalid arguments for GQCP::QCModel::RHF::calculateOrthonormalBasis1DM().
 */
BOOST_AUTO_TEST_CASE(RHF_1DM_invalid_argument) {

    const size_t K = 5;          // The number of spatial orbitals.
    const size_t N_invalid = 3;  // The number of electrons must be even.
    const size_t N_valid = 4;

    BOOST_CHECK_THROW(GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N_invalid), std::invalid_argument);
    BOOST_CHECK_NO_THROW(GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N_valid));
}


/**
 *  Check if GQCP::QCModel::RHF::calculateOrthonormalBasis1DM() yields a correct 1-DM for an example.
 */
BOOST_AUTO_TEST_CASE(RHF_1DM_matrix) {

    const size_t K = 5;  // The number of spatial orbitals.
    const size_t N = 6;  // The number of electrons.
    GQCP::SquareMatrix<double> D_ref {K};
    // clang-format off
    D_ref << 2, 0, 0, 0, 0,
             0, 2, 0, 0, 0,
             0, 0, 2, 0, 0,
             0, 0, 0, 0, 0,
             0, 0, 0, 0, 0;
    // clang-format on

    BOOST_CHECK(GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N).matrix().isApprox(D_ref));
}


/**
 *  Check if the RHF HOMO and LUMO indices are correctly implemented.
 */
BOOST_AUTO_TEST_CASE(HOMO_LUMO_index) {

    // For K=7 and N=10, the index of the HOMO should be 4.
    const size_t K = 7;   // The number of spatial orbitals.
    const size_t N = 10;  // The number of electrons.

    BOOST_CHECK_EQUAL(GQCP::QCModel::RHF<double>::homoIndex(N), 4);
    BOOST_CHECK_EQUAL(GQCP::QCModel::RHF<double>::lumoIndex(K, N), 5);

    BOOST_CHECK_THROW(GQCP::QCModel::RHF<double>::homoIndex(N + 1), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::QCModel::RHF<double>::lumoIndex(K, N + 1), std::invalid_argument);
}


/**
 *  Check if the methods for returning spin-orbital energies are correctly implemented.
 */
BOOST_AUTO_TEST_CASE(spinorbitalEnergies) {

    // Set up toy RHF model parameters.
    const size_t K = 2;
    const auto C = GQCP::RTransformation<double>::Identity(K);
    GQCP::VectorX<double> orbital_energies {K};
    orbital_energies << -0.5, 0.5;

    GQCP::QCModel::RHF<double> rhf_parameters {1, orbital_energies, C};


    // Provide reference values and check the results.
    GQCP::VectorX<double> ref_spinorbital_energies_interleaved {2 * K};
    ref_spinorbital_energies_interleaved << -0.5, -0.5, 0.5, 0.5;
    BOOST_CHECK(rhf_parameters.spinOrbitalEnergiesInterleaved().isApprox(ref_spinorbital_energies_interleaved, 1.0e-12));

    GQCP::VectorX<double> ref_spinorbital_energies_blocked {2 * K};
    ref_spinorbital_energies_blocked << -0.5, 0.5, -0.5, 0.5;
    BOOST_CHECK(rhf_parameters.spinOrbitalEnergiesBlocked().isApprox(ref_spinorbital_energies_blocked, 1.0e-12));
}


/**
 *  Check if the RHF energy is equal to the expectation value of the Hamiltonian through its density matrices.
 */
BOOST_AUTO_TEST_CASE(RHF_DMs) {

    // Perform an RHF calculation.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In an AO basis.

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), hamiltonian, spin_orbital_basis.overlap());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};

    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();
    const auto rhf_energy = rhf_qc_structure.groundStateEnergy();

    // Determine the RHF energy through the expectation value of the Hamiltonian, and check the result.
    // Do the calculations in the RHF MO basis, in order to check the implementation of the RHF density matrices in MO basis.
    hamiltonian.transform(rhf_parameters.expansion());
    const auto D_MO = rhf_parameters.calculateOrthonormalBasis1DM();
    const auto d_MO = rhf_parameters.calculateOrthonormalBasis2DM();
    const double expectation_value = hamiltonian.calculateExpectationValue(D_MO, d_MO);

    BOOST_CHECK(std::abs(rhf_energy - expectation_value) < 1.0e-12);
}


/**
 *  Check the calculation of the ipsocentric current density and the intermediates for its calculation. The test system is H2, 1 au apart in an STO-3G basis set.
 *
 *  The reference implementation is a GAMESS-SYSMO combination. Calculations were performed by Remco Havenith.
 */
BOOST_AUTO_TEST_CASE(ipsocentric_magnetic_inducibility_H2) {

    using namespace GQCP::literals;

    // Set up the molecular Hamiltonian in AO basis.
    const GQCP::Molecule molecule {{GQCP::Nucleus(1, 0.0, 0.0, 0.0), GQCP::Nucleus(1, 0.0, 0.0, 1.0)}};
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

    const auto J_field = GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(grid, orbital_space, x, y, j_op);
    const auto& J_values = J_field.values();

    const auto J_x_field = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/h2_sysmo_currents_x.rgrid");
    const auto& J_x_values = J_x_field.values();

    for (size_t i = 0; i < grid.numberOfPoints(); i++) {
        BOOST_CHECK(J_values[i].col(0).real().isApprox(-J_x_values[i], 1.0e-05));
    }
}


/**
 *  Check the calculation of the ipsocentric current density and the intermediates for its calculation. The test system is H2O in an STO-3G basis set.
 *
 *  The reference implementation is a GAMESS-SYSMO combination. Calculations were performed by Remco Havenith.
 */
BOOST_AUTO_TEST_CASE(ipsocentric_magnetic_inducibility_H2O) {

    using namespace GQCP::literals;

    // Set up the molecular Hamiltonian in AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
    const auto N_P = molecule.numberOfElectronPairs();

    const std::string basis_set {"STO-3G"};
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, basis_set};

    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


    // Read in the GAMESS-UK RHF wave function model parameters. Even though GQCP finds the same orbitals and orbital energies, the phase factors (and orbital coefficients for degenerated orbitals) are unlikely to be reproducible.
    GQCP::MatrixX<double> C_matrix {7, 7};
    // clang-format off
    C_matrix << -0.9944345891,  0.2391588407,  0.0000000000,  0.0936833200,  0.0000000000, -0.1116398665,  0.0000000000,
                -0.0240970451, -0.8857356063,  0.0000000000, -0.4795860800,  0.0000000000,  0.6695790595, -0.0000000000,
                -0.0000000000, -0.0000000000, -0.6072843744, -0.0000000000,  0.0000000000,  0.0000000000,  0.9192346184,
                -0.0031615225, -0.0858966743, -0.0000000000,  0.7474315269,  0.0000000000,  0.7384884432, -0.0000000000,
                 0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  1.0000000000,  0.0000000000,  0.0000000000,
                 0.0045937355, -0.1440394612, -0.4529980986,  0.3294709604,  0.0000000000, -0.7098496636, -0.7324604967,
                 0.0045937355, -0.1440394612,  0.4529980986,  0.3294709604,  0.0000000000, -0.7098496636,  0.7324604967;
    // clang-format on
    GQCP::RTransformation<double> C {C_matrix};

    GQCP::VectorX<double> orbital_energies {7};
    orbital_energies << -20.2628890228, -1.2096967079, -0.5479639306, -0.4365265871, -0.3875856791, 0.4776190290, 0.5881400202;

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

    GQCP::MatrixX<double> A_ref {10, 10};
    // clang-format off
    A_ref << 39.96940,  0.00000,  0.03448,  0.00000,  0.00000, -0.00925,  0.03956,  0.00000,  0.00000,  0.00000,
              0.00000, 40.06171,  0.00000,  0.05605, -0.02686,  0.00000,  0.00000,  0.02998,  0.00000,  0.00000,
              0.03448,  0.00000,  2.37666,  0.00000,  0.00000, -0.11291, -0.02343,  0.00000,  0.00000,  0.00000,
              0.00000,  0.05605,  0.00000,  2.51331, -0.03566,  0.00000,  0.00000, -0.00228,  0.00000,  0.00000,
              0.00000, -0.02686,  0.00000, -0.03566,  1.10153,  0.00000,  0.00000, -0.07674,  0.00000,  0.00000,
             -0.00925,  0.00000, -0.11291,  0.00000,  0.00000,  1.42443,  0.01478,  0.00000,  0.00000,  0.00000,
              0.03956,  0.00000, -0.02343,  0.00000,  0.00000,  0.01478,  0.93180,  0.00000,  0.00000,  0.00000,
              0.00000,  0.02998,  0.00000, -0.00228, -0.07674,  0.00000,  0.00000,  1.03049,  0.00000,  0.00000,
              0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.64372,  0.00000,
              0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.78206;
    // clang-format on

    BOOST_CHECK(A_ref.isApprox((A.asMatrix() * (0.5_ii)).real(), 1.0e-04));


    // Solve the CPHF equations for the angular momentum operator.
    const auto L = complex_spin_orbital_basis.quantize(GQCP::AngularMomentumOperator());
    const auto F_B = rhf_parameters.calculateMagneticFieldResponseForce(L);

    GQCP::MatrixX<double> F_B_ref {3, 10};
    // clang-format off
    F_B_ref << 0.00000, 0.00000, 0.00000,  0.00000, 0.00000,  0.00000, 0.00000, 0.00000,  0.42752, -0.00000,
               0.00000, 0.00000, 0.00000,  0.00000, 0.00000,  0.00000, 0.00000, 0.00000, -0.00000, -0.52599,
               0.00000, 0.12303, 0.00000, -0.06987, 0.26617, -0.00000, 0.00000, 0.47815,  0.00000,  0.00000;
    // clang-format on

    BOOST_CHECK(F_B_ref.transpose().isApprox((F_B * (0.5_ii)).real(), 1.0e-04));


    auto environment_B = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_B);
    auto solver_B = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_B.perform(environment_B);
    const auto x = environment_B.x;

    GQCP::MatrixX<double> x_ref {3, 10};
    // clang-format off
    x_ref << 0.00000, 0.00000, 0.00000,  0.00000, 0.00000,  0.00000, 0.00000, 0.00000,  0.66414, -0.00000,
             0.00000, 0.00000, 0.00000,  0.00000, 0.00000,  0.00000, 0.00000, 0.00000, -0.00000, -0.67257,
             0.00000, 0.00293, 0.00000, -0.02353, 0.27468, -0.00000, 0.00000, 0.48432,  0.00000,  0.00000;
    // clang-format on

    BOOST_CHECK(x_ref.transpose().isApprox(-x.real(), 1.0e-04));


    // Solve the CPHF equations for the linear momentum operator.
    const auto p = complex_spin_orbital_basis.quantize(GQCP::LinearMomentumOperator());
    const auto F_G = rhf_parameters.calculateGaugeOriginTranslationResponseForce(p);

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
    GQCP::MatrixX<GQCP::complex> F_G_reduced {10, 3};
    F_G_reduced.col(0) = 2 * F_G.col(3);  // x <--> yz
    F_G_reduced.col(1) = 2 * F_G.col(4);  // y <--> zx
    F_G_reduced.col(2) = 2 * F_G.col(0);  // z <--> xy

    GQCP::MatrixX<double> F_G_ref {3, 10};
    // clang-format off
    F_G_ref << -0.00000, -1.73412, 0.00000, 0.12306, 0.84094, 0.00000, -0.00000, -0.67978,  0.00000, 0.00000,
               -1.39805,  0.00000, 0.06859, 0.00000, 0.00000, 0.61259, -0.64447,  0.00000,  0.00000, 0.00000,
                0.00000,  0.00000, 0.00000, 0.00000, 0.00000, 0.00000,  0.00000,  0.00000, -0.18459, 0.00000;
    // clang-format on

    BOOST_CHECK(F_G_ref.transpose().isApprox((F_G_reduced * (0.5_ii)).real(), 1.0e-04));


    auto environment_G = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_G);
    auto solver_G = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_G.perform(environment_G);

    const auto y = environment_G.x;

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
    GQCP::MatrixX<GQCP::complex> y_reduced {10, 3};
    y_reduced.col(0) = 2 * y.col(3);  // x <--> yz
    y_reduced.col(1) = 2 * y.col(4);  // y <--> zx
    y_reduced.col(2) = 2 * y.col(0);  // z <--> xy

    GQCP::MatrixX<double> y_ref {3, 10};
    // clang-format off
    y_ref << -0.00000, -0.04243, 0.00000, 0.05961, 0.72221, 0.00000, -0.00000, -0.60452,  0.00000, 0.00000,
             -0.03422,  0.00000, 0.04342, 0.00000, 0.00000, 0.44050, -0.69608,  0.00000,  0.00000, 0.00000,
              0.00000,  0.00000, 0.00000, 0.00000, 0.00000, 0.00000,  0.00000,  0.00000, -0.28675, 0.00000;
    // clang-format on

    BOOST_CHECK(y_ref.transpose().isApprox(-y_reduced.real(), 1.0e-04));


    // Calculate the ipsocentric magnetic inducibility on a grid and check the results.
    const auto j_op = complex_spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());

    const GQCP::Vector<double, 3> origin {-2.0, -2.0, -2.0};
    const std::array<size_t, 3> steps {6, 6, 6};
    const std::array<double, 3> step_sizes {0.8, 0.8, 0.8};
    const GQCP::CubicGrid grid {origin, steps, step_sizes};

    const auto J_field = GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(grid, orbital_space, x, y, j_op);
    const auto& J_values = J_field.values();

    const auto J_x_field = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/h2o_sysmo_currents_x.rgrid");
    const auto& J_x_values = J_x_field.values();

    const auto J_y_field = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/h2o_sysmo_currents_y.rgrid");
    const auto& J_y_values = J_y_field.values();

    const auto J_z_field = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/h2o_sysmo_currents_z.rgrid");
    const auto& J_z_values = J_z_field.values();

    for (size_t i = 0; i < grid.numberOfPoints(); i++) {
        BOOST_CHECK(J_values[i].col(0).real().isApprox(-J_x_values[i], 1.0e-05));
        BOOST_CHECK(J_values[i].col(1).real().isApprox(-J_y_values[i], 1.0e-05));
        BOOST_CHECK(J_values[i].col(2).real().isApprox(-J_z_values[i], 1.0e-05));
    }
}


/**
 *  Check the calculation of the ipsocentric current density and the intermediates for its calculation. The test system is CH4 in an STO-3G basis set.
 *
 *  The reference implementation is a GAMESS-SYSMO combination. Calculations were performed by Remco Havenith.
 */
BOOST_AUTO_TEST_CASE(ipsocentric_magnetic_inducibility_CH4) {

    using namespace GQCP::literals;

    // Set up the molecular Hamiltonian in AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/ch4_crawdad.xyz");
    const auto N_P = molecule.numberOfElectronPairs();

    const std::string basis_set {"STO-3G"};
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, basis_set};

    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


    // Read in the GAMESS-UK RHF wave function model parameters. Even though GQCP finds the same orbitals and orbital energies, the phase factors (and orbital coefficients for degenerated orbitals) are unlikely to be reproducible.
    GQCP::MatrixX<double> C_matrix {9, 9};
    // clang-format off
    C_matrix <<  0.9919348327,  0.2213540576, -0.0000000000,  0.0000000000, -0.0000000000, -0.0000000000,  0.0000000000,  0.0000000000,  0.2521385875,
                 0.0382620669, -0.6281391692, -0.0000000000,  0.0000000000, -0.0000000000, -0.0000000000,  0.0000000000,  0.0000000000, -1.6275938477,
                 0.0000000000, -0.0000000000, -0.0170203876,  0.5712733267,  0.0149068523,  0.9277617097, -0.0501500162,  0.5980861308, -0.0000000000,
                -0.0000000000,  0.0000000000, -0.0362386757, -0.0159623991,  0.5703482109, -0.4350188268, -0.8147978414,  0.6064873718,  0.0000000000,
                 0.0000000000, -0.0000000000,  0.5703176193,  0.0160346393,  0.0366854952,  0.4134980825, -0.7446833471, -0.7038678014, -0.0000000000,
                -0.0069828775, -0.1806242559, -0.2901680552,  0.3007521493, -0.3117701078, -0.5270666072, -0.8380203988, -0.3861413420,  0.6663725071,
                -0.0069828775, -0.1806242559,  0.2722447054,  0.3008282225,  0.3274677947, -0.5031688400,  0.8937096170, -0.2780049862,  0.6663725071,
                -0.0069828775, -0.1806242559, -0.3104059065, -0.3176374550,  0.2731383827,  0.9862366485,  0.0110848034, -0.3954705220,  0.6663725071,
                -0.0069828775, -0.1806242559,  0.3283292563, -0.2839429169, -0.2888360696,  0.0439987988, -0.0667740216,  1.0596168501,  0.6663725071;
    // clang-format on
    GQCP::RTransformation<double> C {C_matrix};

    GQCP::VectorX<double> orbital_energies {9};
    orbital_energies << -11.0298571295, -0.9110637967, -0.5197078520, -0.5197078519, -0.5197078518, 0.7174507706, 0.717450770, 0.7174507706, 0.7580377307;

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

    GQCP::MatrixX<double> A_ref {20, 20};
    // clang-format off
    A_ref << 22.0723066626,  0.0000000000,  0.0000000000,  0.0000000000, -0.0454611397,  0.0000000000,  0.0000000000,  0.0000000000,  0.0069922537,  0.0056960124, -0.0025782743,  0.0145049104,  0.0033466260,  0.0001134245, -0.0047798183,  0.0334377141, -0.0058404038,  0.0062809647,  0.0030751869, -0.0134784855,
              0.0000000000, 22.0723066626,  0.0000000000,  0.0000000000,  0.0000000000, -0.0454611397,  0.0000000000,  0.0000000000,  0.0056960124, -0.0003372180,  0.0040976683, -0.0242564874,  0.0001134245, -0.0100579300, -0.0009844264, -0.0016967976,  0.0062809647, -0.0009238452,  0.0036419263, -0.0303131393,
              0.0000000000,  0.0000000000, 22.0723066626,  0.0000000000,  0.0000000000,  0.0000000000, -0.0454611397,  0.0000000000, -0.0025782743,  0.0040976683, -0.0066550357, -0.0266715787, -0.0047798183, -0.0009844264,  0.0067113040,  0.0197277185,  0.0030751869,  0.0036419263,  0.0067642491,  0.0202382492,
              0.0000000000,  0.0000000000,  0.0000000000, 22.2839653047,  0.0000000000,  0.0000000000,  0.0000000000, -0.0026696604,  0.0031495210, -0.0052669278, -0.0057913282,  0.0000000000,  0.0072604918, -0.0003684338,  0.0042835744,  0.0000000000, -0.0029266484, -0.0065820378,  0.0043944284,  0.0000000000,
             -0.0454611397,  0.0000000000,  0.0000000000,  0.0000000000,  2.3644122877,  0.0000000000,  0.0000000000,  0.0000000000, -0.0377539348, -0.0307550170,  0.0139211195, -0.0187461840, -0.0180697536, -0.0006124238,  0.0258081238, -0.0432149887,  0.0315346433, -0.0339134052, -0.0166041465,  0.0174196297,
              0.0000000000, -0.0454611397,  0.0000000000,  0.0000000000,  0.0000000000,  2.3644122878,  0.0000000000,  0.0000000000, -0.0307550170,  0.0018207729, -0.0221249270,  0.0313491475, -0.0006124237,  0.0543067302,  0.0053153063,  0.0021929456, -0.0339134052,  0.0049882047, -0.0196641959,  0.0391767802,
              0.0000000000,  0.0000000000, -0.0454611397,  0.0000000000,  0.0000000000,  0.0000000000,  2.3644122878,  0.0000000000,  0.0139211195, -0.0221249271,  0.0359331620,  0.0344704178,  0.0258081238,  0.0053153063, -0.0362369766, -0.0254961547, -0.0166041465, -0.0196641959, -0.0365228480, -0.0261559660,
              0.0000000000,  0.0000000000,  0.0000000000, -0.0026696604,  0.0000000000,  0.0000000000,  0.0000000000,  2.5429918963, -0.0412439043,  0.0689719700,  0.0758391476,  0.0000000000, -0.0950782763,  0.0048247493, -0.0560946679,  0.0000000000,  0.0383253223,  0.0861937222, -0.0575463337,  0.0000000000,
              0.0069922537,  0.0056960124, -0.0025782743,  0.0031495210, -0.0377539348, -0.0307550170,  0.0139211195, -0.0412439043,  1.6257891339, -0.0098766123, -0.0102622786, -0.0049745885,  0.0077855755, -0.0175001064, -0.0210476922, -0.0269171691, -0.0032698099,  0.0087151220,  0.0062703511,  0.0495679722,
              0.0056960124, -0.0003372180,  0.0040976683, -0.0052669278, -0.0307550170,  0.0018207729, -0.0221249271,  0.0689719700, -0.0098766123,  1.6371286000,  0.0176149098,  0.0028724975,  0.0034440377,  0.0009033101,  0.0035699364, -0.0422331308, -0.0117006844,  0.0121279745,  0.0207966036, -0.0025442187,
             -0.0025782743,  0.0040976683, -0.0066550357, -0.0057913282,  0.0139211195, -0.0221249270,  0.0359331620,  0.0758391476, -0.0102622786,  0.0176149098,  1.6383823499, -0.0063236158,  0.0103227465, -0.0103786685, -0.0086888856,  0.0348670663,  0.0045143705, -0.0138073926, -0.0088581646,  0.0348196942,
              0.0145049104, -0.0242564874, -0.0266715787,  0.0000000000, -0.0187461840,  0.0313491475,  0.0344704178,  0.0000000000, -0.0049745885,  0.0028724975, -0.0063236158,  1.7274849107, -0.0269171691, -0.0422331308,  0.0348670663,  0.0000000000,  0.0495679722, -0.0025442186,  0.0348196942,  0.0000000000,
              0.0033466260,  0.0001134245, -0.0047798183,  0.0072604918, -0.0180697536, -0.0006124237,  0.0258081238, -0.0950782763,  0.0077855755,  0.0034440377,  0.0103227465, -0.0269171691,  1.6499615391, -0.0015560891,  0.0182751951, -0.0027054652, -0.0072554515,  0.0054186840, -0.0088579260,  0.0225682579,
              0.0001134245, -0.0100579300, -0.0009844264, -0.0003684338, -0.0006124238,  0.0543067302,  0.0053153063,  0.0048247493, -0.0175001064,  0.0009033101, -0.0103786685, -0.0422331308, -0.0015560891,  1.6201186871, -0.0009515486, -0.0002234507, -0.0221831792,  0.0006239925, -0.0125732182, -0.0440295251,
             -0.0047798183, -0.0009844264,  0.0067113040,  0.0042835744,  0.0258081238,  0.0053153063, -0.0362369766, -0.0560946679, -0.0210476922,  0.0035699364, -0.0086888856,  0.0348670663,  0.0182751951, -0.0009515486,  1.6312198571,  0.0040011679,  0.0162446094,  0.0024376126,  0.0066314590, -0.0357428279,
              0.0334377141, -0.0016967976,  0.0197277185,  0.0000000000, -0.0432149887,  0.0021929456, -0.0254961547,  0.0000000000, -0.0269171691, -0.0422331308,  0.0348670663,  0.0000000000, -0.0027054652, -0.0002234507,  0.0040011679,  1.7274849105,  0.0225682579, -0.0440295251, -0.0357428279,  0.0000000000,
             -0.0058404038,  0.0062809647,  0.0030751869, -0.0029266484,  0.0315346433, -0.0339134052, -0.0166041465,  0.0383253223, -0.0032698099, -0.0117006844,  0.0045143705,  0.0495679722, -0.0072554515, -0.0221831792,  0.0162446094,  0.0225682579,  1.6255494103,  0.0114327012, -0.0080129165,  0.0076800536,
              0.0062809647, -0.0009238452,  0.0036419263, -0.0065820378, -0.0339134052,  0.0049882047, -0.0196641959,  0.0861937222,  0.0087151220,  0.0121279745, -0.0138073926, -0.0025442186,  0.0054186840,  0.0006239925,  0.0024376126, -0.0440295251,  0.0114327012,  1.6440527963, -0.0166633610, -0.0026490468,
              0.0030751869,  0.0036419263,  0.0067642491,  0.0043944284, -0.0166041465, -0.0196641959, -0.0365228480, -0.0575463337,  0.0062703511,  0.0207966036, -0.0088581646,  0.0348196942, -0.0088579260, -0.0125732182,  0.0066314590, -0.0357428279, -0.0080129165, -0.0166633610,  1.6316978765,  0.0023224479,
             -0.0134784855, -0.0303131393,  0.0202382492,  0.0000000000,  0.0174196297,  0.0391767802, -0.0261559660,  0.0000000000,  0.0495679722, -0.0025442187,  0.0348196942,  0.0000000000,  0.0225682579, -0.0440295251, -0.0357428279,  0.0000000000,  0.0076800536, -0.0026490468,  0.0023224479,  1.7274849104;
    // clang-format on

    BOOST_CHECK(A_ref.isApprox((A.asMatrix() * (0.5_ii)).real(), 1.0e-06));


    // Solve the CPHF equations for the angular momentum operator.
    const auto L = complex_spin_orbital_basis.quantize(GQCP::AngularMomentumOperator());
    const auto F_B = rhf_parameters.calculateMagneticFieldResponseForce(L);

    GQCP::MatrixX<double> F_B_ref {3, 20};
    // clang-format off
    F_B_ref << -0.00000, -0.00000,  0.00000,  0.00000, -0.00000, -0.00000,  0.00000,  0.00000, -0.16201, -0.34170,  0.22266, -0.00000, -0.00026, -0.01734, -0.00105, -0.00000, -0.17499,  0.27440,  0.29446,  0.00000,
               -0.00000,  0.00000, -0.00000, -0.00000, -0.00000,  0.00000,  0.00000,  0.00000, -0.37261,  0.02869, -0.22873, -0.00000,  0.15383, -0.29509, -0.28611, -0.00000, -0.01937, -0.00644, -0.02254, -0.00000,
               -0.00000, -0.00000,  0.00000, -0.00000,  0.00000, -0.00000,  0.00000, -0.00000, -0.02851, -0.00837, -0.00789, -0.00000,  0.16242,  0.32404, -0.24742,  0.00000,  0.37225, -0.01144,  0.23078,  0.00000;
    // clang-format on

    BOOST_CHECK(F_B_ref.transpose().isApprox((F_B * (0.5_ii)).real(), 1.0e-04));


    auto environment_B = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_B);
    auto solver_B = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_B.perform(environment_B);
    const auto x = environment_B.x;

    GQCP::MatrixX<double> x_ref {3, 20};
    // clang-format off
    x_ref << 0.00000, -0.00000,  0.00000, -0.00000, -0.00000,  0.00000, -0.00000,  0.00000, -0.10206, -0.21526,  0.14027, -0.00000, -0.00016, -0.01092, -0.00066, -0.00000, -0.11024,  0.17286,  0.18550, -0.00000,
             0.00000,  0.00000, -0.00000,  0.00000, -0.00000, -0.00000,  0.00000, -0.00000, -0.23473,  0.01807, -0.14409,  0.00000,  0.09691, -0.18590, -0.18024, -0.00000, -0.01220, -0.00405, -0.01420,  0.00000,
             0.00000,  0.00000,  0.00000,  0.00000, -0.00000, -0.00000, -0.00000, -0.00000, -0.01796, -0.00528, -0.00497,  0.00000,  0.10232,  0.20414, -0.15587, -0.00000,  0.23450, -0.00720,  0.14539,  0.00000;
    // clang-format on

    BOOST_CHECK(x_ref.transpose().isApprox(-x.real(), 1.0e-04));


    // Solve the CPHF equations for the linear momentum operator.
    const auto p = complex_spin_orbital_basis.quantize(GQCP::LinearMomentumOperator());
    const auto F_G = rhf_parameters.calculateGaugeOriginTranslationResponseForce(p);

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
    GQCP::MatrixX<GQCP::complex> F_G_reduced {20, 3};
    F_G_reduced.col(0) = 2 * F_G.col(3);  // x <--> yz
    F_G_reduced.col(1) = 2 * F_G.col(4);  // y <--> zx
    F_G_reduced.col(2) = 2 * F_G.col(0);  // z <--> xy

    GQCP::MatrixX<double> F_G_ref {3, 20};
    // clang-format off
    F_G_ref <<  1.26367, -0.06831,  0.81463, -0.00000,  0.17853, -0.00965,  0.11509, -0.00000,  0.31054,  0.51667, -0.43840, -0.02284,  0.01602,  0.00139, -0.02474,  0.76660, -0.25955, 0.53664,  0.44761, 0.02000,
               -0.59252, -1.10981,  0.82608,  0.00000, -0.08371, -0.15679,  0.11671,  0.00000, -0.61627,  0.01880, -0.41678, -0.04863, -0.29640,  0.50311,  0.46332, -0.02142, -0.04745, 0.01528, -0.01351, 0.76536,
                0.56321, -1.01431, -0.95871, -0.00000,  0.07957, -0.14330, -0.13545,  0.00000,  0.03095, -0.01852,  0.03777,  0.76532,  0.31083,  0.54850, -0.39771,  0.02152, -0.61695, 0.04810, -0.41333, 0.04923;
    // clang-format on

    BOOST_CHECK(F_G_ref.transpose().isApprox((F_G_reduced * (0.5_ii)).real(), 1.0e-04));


    auto environment_G = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_G);
    auto solver_G = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_G.perform(environment_G);

    const auto y = environment_G.x;

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
    GQCP::MatrixX<GQCP::complex> y_reduced {20, 3};
    y_reduced.col(0) = 2 * y.col(3);  // x <--> yz
    y_reduced.col(1) = 2 * y.col(4);  // y <--> zx
    y_reduced.col(2) = 2 * y.col(0);  // z <--> xy

    GQCP::MatrixX<double> y_ref {3, 20};
    // clang-format off
    y_ref << 0.05639, -0.00305,  0.03635,  0.00000,  0.10350, -0.00559,  0.06672, -0.00000,  0.19682,  0.32746, -0.27785, -0.01428,  0.01016, 0.00088, -0.01568,  0.47926, -0.16450, 0.34011,  0.28369, 0.01251,
            -0.02644, -0.04953,  0.03686,  0.00000, -0.04853, -0.09090,  0.06766, -0.00000, -0.39058,  0.01192, -0.26415, -0.03040, -0.18785, 0.31887,  0.29365, -0.01339, -0.03007, 0.00968, -0.00856, 0.47848,
             0.02513, -0.04526, -0.04278, -0.00000,  0.04613, -0.08307, -0.07852,  0.00000,  0.01961, -0.01173,  0.02394,  0.47846,  0.19700, 0.34763, -0.25206,  0.01345, -0.39102, 0.03049, -0.26196, 0.03078;
    // clang-format on

    BOOST_CHECK(y_ref.transpose().isApprox(-y_reduced.real(), 1.0e-04));


    // Calculate the ipsocentric magnetic inducibility on a grid and check the results.
    const auto j_op = complex_spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());

    const GQCP::Vector<double, 3> origin {-2.0, -2.0, -2.0};
    const std::array<size_t, 3> steps {6, 6, 6};
    const std::array<double, 3> step_sizes {0.8, 0.8, 0.8};

    const GQCP::CubicGrid grid {origin, steps, step_sizes};

    const auto J_field = GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(grid, orbital_space, x, y, j_op);
    const auto& J_values = J_field.values();

    const auto J_x_field = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/ch4_sysmo_currents_x.rgrid");
    const auto& J_x_values = J_x_field.values();

    const auto J_y_field = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/ch4_sysmo_currents_y.rgrid");
    const auto& J_y_values = J_y_field.values();

    const auto J_z_field = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/ch4_sysmo_currents_z.rgrid");
    const auto& J_z_values = J_z_field.values();

    for (size_t i = 0; i < grid.numberOfPoints(); i++) {
        BOOST_CHECK(J_values[i].col(0).real().isApprox(-J_x_values[i], 1.0e-06));
        BOOST_CHECK(J_values[i].col(1).real().isApprox(-J_y_values[i], 1.0e-06));
        BOOST_CHECK(J_values[i].col(2).real().isApprox(-J_z_values[i], 1.0e-06));
    }
}


/**
 *  Check the calculation of the ipsocentric current density and the intermediates for its calculation. The test system is CO in an STO-3G basis set.
 *
 *  The reference implementation is a GAMESS-SYSMO combination. Calculations were performed by Remco Havenith.
 */
BOOST_AUTO_TEST_CASE(ipsocentric_magnetic_inducibility_CO) {

    using namespace GQCP::literals;

    // Set up the molecular Hamiltonian in AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/CO_mulliken.xyz");
    const auto N_P = molecule.numberOfElectronPairs();

    const std::string basis_set {"STO-3G"};
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, basis_set};

    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


    // Read in the GAMESS-UK RHF wave function model parameters. Even though GQCP finds the same orbitals and orbital energies, the phase factors (and orbital coefficients for degenerated orbitals) are unlikely to be reproducible.
    GQCP::MatrixX<double> C_matrix {10, 10};
    // clang-format off
    C_matrix << -0.0004232010, -0.9936263296, -0.1238367752, -0.1695802599,  0.0000000000,  0.0000000000, -0.1650756758,  0.0000000000,  0.0000000000, -0.1224277929,
                 0.0083068645, -0.0262032651,  0.2436766532,  0.5589623932,  0.0000000000,  0.0000000000,  0.7476988939,  0.0000000000,  0.0000000000,  0.9373564955,
                 0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000, -0.4456287928,  0.0000000000,  0.0000000000, -0.9290557653,  0.0000000000,
                 0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000, -0.4456287928,  0.0000000000,  0.0000000000, -0.9290557653,  0.0000000000,  0.0000000000,
                 0.0071447865, -0.0068485366,  0.1658881312,  0.0648680512,  0.0000000000,  0.0000000000, -0.5747190050,  0.0000000000,  0.0000000000,  1.2058266243,
                -0.9941841695,  0.0001292066, -0.2225383712,  0.1316983515,  0.0000000000,  0.0000000000, -0.0014152641,  0.0000000000,  0.0000000000,  0.1264210162,
                -0.0273361012,  0.0068489345,  0.7705814626, -0.6425157474,  0.0000000000,  0.0000000000,  0.0491584905,  0.0000000000,  0.0000000000, -1.0410108471,
                 0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000, -0.7941917491,  0.0000000000,  0.0000000000,  0.6564976022,  0.0000000000,
                 0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000, -0.7941917491,  0.0000000000,  0.0000000000,  0.6564976022,  0.0000000000,  0.0000000000,
                 0.0065508084, -0.0012098766, -0.2104812688, -0.6146979958,  0.0000000000,  0.0000000000,  0.4445924885,  0.0000000000,  0.0000000000,  0.9584546795;
    // clang-format on
    GQCP::RTransformation<double> C {C_matrix};

    GQCP::VectorX<double> orbital_energies {10};
    orbital_energies << -20.4155574135, -11.0921863361, -1.4452522251, -0.6968296505, -0.5399124938, -0.5399124938, -0.4451302996, 0.3061360876, 0.3061360876, 1.0090710424;

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

    GQCP::MatrixX<double> A_ref {21, 21};
    // clang-format off
    A_ref << 40.1888500628,  0.0000000000,  0.0000000000, -0.0006095647,  0.0000000000,  0.0000000000, -0.0256630925,  0.0000000000,  0.0000000000,  0.0267699982,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0293492009,  0.0000000000,  0.0000000000,  0.0000000000, -0.0042034850,  0.0000000000,  0.0000000000,
              0.0000000000, 40.1888500628,  0.0000000000,  0.0000000000, -0.0006095647,  0.0000000000,  0.0000000000, -0.0256630925,  0.0000000000,  0.0000000000,  0.0267699982,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0293492009,  0.0000000000, -0.0042034850,  0.0000000000,
              0.0000000000,  0.0000000000, 41.1790576544,  0.0000000000,  0.0000000000, -0.0008823284,  0.0000000000,  0.0000000000, -0.0264582361,  0.0000000000,  0.0000000000,  0.0531723303,  0.0027065514,  0.0000000000,  0.0000000000,  0.0000000000,  0.0027065514,  0.0000000000,  0.0000000000,  0.0000000000, -0.0203976501,
             -0.0006095647,  0.0000000000,  0.0000000000, 21.3325791480,  0.0000000000,  0.0000000000, -0.0229280832,  0.0000000000,  0.0000000000, -0.0378886021,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0240231629,  0.0000000000,  0.0000000000,  0.0000000000, -0.0425838536,  0.0000000000,  0.0000000000,
              0.0000000000, -0.0006095647,  0.0000000000,  0.0000000000, 21.3325791480,  0.0000000000,  0.0000000000, -0.0229280832,  0.0000000000,  0.0000000000, -0.0378886021,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0240231629,  0.0000000000, -0.0425838536,  0.0000000000,
              0.0000000000,  0.0000000000, -0.0008823284,  0.0000000000,  0.0000000000, 22.6213676532,  0.0000000000,  0.0000000000, -0.0055434630,  0.0000000000,  0.0000000000, -0.0268446904,  0.0007553989,  0.0000000000,  0.0000000000,  0.0000000000,  0.0007553989,  0.0000000000,  0.0000000000,  0.0000000000, -0.0632430742,
             -0.0256630925,  0.0000000000,  0.0000000000, -0.0229280832,  0.0000000000,  0.0000000000,  2.4698470661,  0.0000000000,  0.0000000000, -0.0900715815,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0672149516,  0.0000000000,  0.0000000000,  0.0000000000,  0.0746442687,  0.0000000000,  0.0000000000,
              0.0000000000, -0.0256630925,  0.0000000000,  0.0000000000, -0.0229280832,  0.0000000000,  0.0000000000,  2.4698470661,  0.0000000000,  0.0000000000, -0.0900715815,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0672149516,  0.0000000000,  0.0746442687,  0.0000000000,
              0.0000000000,  0.0000000000, -0.0264582361,  0.0000000000,  0.0000000000, -0.0055434630,  0.0000000000,  0.0000000000,  3.7301232990,  0.0000000000,  0.0000000000, -0.0528300176, -0.0513834390,  0.0000000000,  0.0000000000,  0.0000000000, -0.0513834390,  0.0000000000,  0.0000000000,  0.0000000000,  0.0333320385,
              0.0267699982,  0.0000000000,  0.0000000000, -0.0378886021,  0.0000000000,  0.0000000000, -0.0900715815,  0.0000000000,  0.0000000000,  1.0979326390,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0929779832,  0.0000000000,  0.0000000000,  0.0000000000,  0.0348165921,  0.0000000000,  0.0000000000,
              0.0000000000,  0.0267699982,  0.0000000000,  0.0000000000, -0.0378886021,  0.0000000000,  0.0000000000, -0.0900715815,  0.0000000000,  0.0000000000,  1.0979326390,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0929779832,  0.0000000000,  0.0348165921,  0.0000000000,
              0.0000000000,  0.0000000000,  0.0531723303,  0.0000000000,  0.0000000000, -0.0268446904,  0.0000000000,  0.0000000000, -0.0528300176,  0.0000000000,  0.0000000000,  2.3812958673, -0.0070214587,  0.0000000000,  0.0000000000,  0.0000000000, -0.0070214587,  0.0000000000,  0.0000000000,  0.0000000000,  0.0239567932,
              0.0000000000,  0.0000000000,  0.0027065514,  0.0000000000,  0.0000000000,  0.0007553989,  0.0000000000,  0.0000000000, -0.0513834390,  0.0000000000,  0.0000000000, -0.0070214587,  0.9151801006,  0.0000000000,  0.0000000000,  0.0000000000, -0.0069219852,  0.0000000000,  0.0000000000,  0.0000000000,  0.0083425924,
              0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.6978343079,  0.0000000000,  0.2242677779,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,
              0.0293492009,  0.0000000000,  0.0000000000,  0.0240231629,  0.0000000000,  0.0000000000,  0.0672149516,  0.0000000000,  0.0000000000,  0.0929779832,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  1.9479123926,  0.0000000000,  0.0000000000,  0.0000000000, -0.1468870360,  0.0000000000,  0.0000000000,
              0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.2242677779,  0.0000000000,  0.6978343079,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,
              0.0000000000,  0.0000000000,  0.0027065514,  0.0000000000,  0.0000000000,  0.0007553989,  0.0000000000,  0.0000000000, -0.0513834390,  0.0000000000,  0.0000000000, -0.0070214587, -0.0069219852,  0.0000000000,  0.0000000000,  0.0000000000,  0.9151801006,  0.0000000000,  0.0000000000,  0.0000000000,  0.0083425924,
              0.0000000000,  0.0293492009,  0.0000000000,  0.0000000000,  0.0240231629,  0.0000000000,  0.0000000000,  0.0672149516,  0.0000000000,  0.0000000000,  0.0929779832,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  1.9479123926,  0.0000000000, -0.1468870360,  0.0000000000,
             -0.0042034850,  0.0000000000,  0.0000000000, -0.0425838536,  0.0000000000,  0.0000000000,  0.0746442687,  0.0000000000,  0.0000000000,  0.0348165921,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000, -0.1468870360,  0.0000000000,  0.0000000000,  0.0000000000,  0.5481901617,  0.0000000000,  0.0000000000,
              0.0000000000, -0.0042034850,  0.0000000000,  0.0000000000, -0.0425838536,  0.0000000000,  0.0000000000,  0.0746442687,  0.0000000000,  0.0000000000,  0.0348165921,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000, -0.1468870360,  0.0000000000,  0.5481901617,  0.0000000000,
              0.0000000000,  0.0000000000, -0.0203976501,  0.0000000000,  0.0000000000, -0.0632430742,  0.0000000000,  0.0000000000,  0.0333320385,  0.0000000000,  0.0000000000,  0.0239567932,  0.0083425924,  0.0000000000,  0.0000000000,  0.0000000000,  0.0083425924,  0.0000000000,  0.0000000000,  0.0000000000,  1.9740177310;
    // clang-format on

    BOOST_CHECK(A_ref.isApprox((A.asMatrix() * (0.5_ii)).real(), 1.0e-06));


    // Solve the CPHF equations for the angular momentum operator.
    const auto L = complex_spin_orbital_basis.quantize(GQCP::AngularMomentumOperator());
    const auto F_B = rhf_parameters.calculateMagneticFieldResponseForce(L);

    GQCP::MatrixX<double> F_B_ref {3, 21};
    // clang-format off
    F_B_ref << -0.58298,  0.00000,  0.00000, -0.86912,  0.00000,  0.00000, -0.26760,  0.00000,  0.00000, -0.42089,  0.00000,  0.00000,  0.00000,  0.00000,  0.38323,  0.00000,  0.00000,  0.00000,  1.05047,  0.00000,  0.00000,
                0.00000,  0.58298,  0.00000,  0.00000,  0.86912,  0.00000,  0.00000,  0.26760,  0.00000,  0.00000,  0.42089,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000, -0.38323,  0.00000, -1.05047,  0.00000,
                0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000;
    // clang-format on

    BOOST_CHECK(F_B_ref.transpose().isApprox((F_B * (0.5_ii)).real(), 1.0e-04));


    auto environment_B = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_B);
    auto solver_B = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_B.perform(environment_B);
    const auto x = environment_B.x;

    GQCP::MatrixX<double> x_ref {3, 21};
    // clang-format off
    x_ref << -0.01437,  0.00000,  0.00000, -0.03813,  0.00000,  0.00000, -0.20024,  0.00000,  0.00000, -0.49912,  0.00000,  0.00000,  0.00000,  0.00000,  0.38464,  0.00000,  0.00000,  0.00000,  2.07521,  0.00000,  0.00000,
              0.00000,  0.01437,  0.00000,  0.00000,  0.03813,  0.00000,  0.00000,  0.20024,  0.00000,  0.00000,  0.49912,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000, -0.38464,  0.00000, -2.07521,  0.00000,
              0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000;
    // clang-format on

    BOOST_CHECK(x_ref.transpose().isApprox(-x.real(), 1.0e-04));


    // Solve the CPHF equations for the linear momentum operator.
    const auto p = complex_spin_orbital_basis.quantize(GQCP::LinearMomentumOperator());
    const auto F_G = rhf_parameters.calculateGaugeOriginTranslationResponseForce(p);

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
    GQCP::MatrixX<GQCP::complex> F_G_reduced {21, 3};
    F_G_reduced.col(0) = 2 * F_G.col(3);  // x <--> yz
    F_G_reduced.col(1) = 2 * F_G.col(4);  // y <--> zx
    F_G_reduced.col(2) = 2 * F_G.col(0);  // z <--> xy

    GQCP::MatrixX<double> F_G_ref {3, 21};
    // clang-format off
    F_G_ref << 0.00000, -1.25807,  0.00000,  0.00000,  1.41754,  0.00000,  0.00000,  0.02240,  0.00000,  0.00000, -0.55694,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000, -0.24625,  0.00000, -0.26685,  0.00000,
              -1.25807,  0.00000,  0.00000,  1.41754,  0.00000,  0.00000,  0.02240,  0.00000,  0.00000, -0.55694,  0.00000,  0.00000,  0.00000,  0.00000, -0.24625,  0.00000,  0.00000,  0.00000, -0.26685,  0.00000,  0.00000,
               0.00000,  0.00000, -1.59193,  0.00000,  0.00000, -1.48694,  0.00000,  0.00000, -0.84885,  0.00000,  0.00000, -0.80094, -0.52589,  0.00000,  0.00000,  0.00000, -0.52589,  0.00000,  0.00000,  0.00000,  1.17267;
    // clang-format on

    BOOST_CHECK(F_G_ref.transpose().isApprox((F_G_reduced * (0.5_ii)).real(), 1.0e-04));


    auto environment_G = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_G);
    auto solver_G = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_G.perform(environment_G);

    const auto y = environment_G.x;

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
    GQCP::MatrixX<GQCP::complex> y_reduced {21, 3};
    y_reduced.col(0) = 2 * y.col(3);  // x <--> yz
    y_reduced.col(1) = 2 * y.col(4);  // y <--> zx
    y_reduced.col(2) = 2 * y.col(0);  // z <--> xy

    GQCP::MatrixX<double> y_ref {3, 21};
    // clang-format off
    y_ref << 0.00000, -0.03093,  0.00000,  0.00000,  0.06479,  0.00000,  0.00000,  0.01069,  0.00000,  0.00000, -0.47584,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000, -0.14144,  0.00000, -0.49112,  0.00000,
            -0.03093,  0.00000,  0.00000,  0.06479,  0.00000,  0.00000,  0.01069,  0.00000,  0.00000, -0.47584,  0.00000,  0.00000,  0.00000,  0.00000, -0.14144,  0.00000,  0.00000,  0.00000, -0.49112,  0.00000,  0.00000,
             0.00000,  0.00000, -0.03799,  0.00000,  0.00000, -0.06448,  0.00000,  0.00000, -0.25489,  0.00000,  0.00000, -0.35152, -0.60154,  0.00000,  0.00000,  0.00000, -0.60154,  0.00000,  0.00000,  0.00000,  0.60525;
    // clang-format on

    BOOST_CHECK(y_ref.transpose().isApprox(-y_reduced.real(), 1.0e-04));


    // Calculate the ipsocentric magnetic inducibility on a grid and check the results.
    const auto j_op = complex_spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());

    const GQCP::Vector<double, 3> origin {-2.0, -2.0, -2.0};
    const std::array<size_t, 3> steps {6, 6, 6};
    const std::array<double, 3> step_sizes {0.8, 0.8, 0.8};

    const GQCP::CubicGrid grid {origin, steps, step_sizes};

    const auto J_field = GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(grid, orbital_space, x, y, j_op);
    const auto& J_values = J_field.values();

    const auto J_x_field = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/co_sysmo_currents_x.rgrid");
    const auto& J_x_values = J_x_field.values();

    const auto J_z_field = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/co_sysmo_currents_z.rgrid");
    const auto& J_z_values = J_z_field.values();

    for (size_t i = 0; i < grid.numberOfPoints(); i++) {
        BOOST_CHECK(J_values[i].col(0).real().isApprox(-J_x_values[i], 1.0e-06));
        BOOST_CHECK(J_values[i].col(2).real().isApprox(-J_z_values[i], 1.0e-04));
    }
}
