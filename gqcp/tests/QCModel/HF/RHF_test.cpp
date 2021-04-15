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
BOOST_AUTO_TEST_CASE(ipsocentric_H2) {

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
    const auto grid_points = grid.points();

    const auto J_field = GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(grid, orbital_space, x, y, j_op);
}


/**
 *  Check the calculation of the ipsocentric current density and the intermediates for its calculation. The test system is H2O in an STO-3G basis set.
 *
 *  The reference implementation is a GAMESS-SYSMO combination. Calculations were performed by Remco Havenith.
 */
BOOST_AUTO_TEST_CASE(ipsocentric_H2O) {

    using namespace GQCP::literals;

    // Set up the molecular Hamiltonian in AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
    const auto N_P = molecule.numberOfElectronPairs();

    const std::string basis_set {"STO-3G"};
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, basis_set};

    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


    // Read in the GAMESS-UK RHF wave function model parameters. Even though GQCP finds the same orbitals and orbital energies, the phase factors (and orbital coefficients for degenerated orbitals) are unlikely to be reproducible.
    GQCP::MatrixX<double> C_matrix {2, 2};
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
    GQCP::RTransformation<GQCP::complex> C_adjusted_complex {C_adjusted.matrix().cast<GQCP::complex>()};
    complex_spin_orbital_basis.transform(C_adjusted_complex);

    spin_orbital_basis.transform(C_adjusted);
    hamiltonian.transform(C_adjusted);


    // Calculate the orbital Hessian.
    const auto orbital_space = rhf_parameters.orbitalSpace();
    auto A = rhf_parameters_adjusted.calculateOrbitalHessianForImaginaryResponse(hamiltonian, orbital_space);

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
    const auto F_B = rhf_parameters_adjusted.calculateMagneticFieldResponseForce(L);

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
    const auto F_G = rhf_parameters_adjusted.calculateGaugeOriginTranslationResponseForce(p);

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
    // std::cout << y << std::endl;

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
}
