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

#define BOOST_TEST_MODULE "RHFSCFSolver"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"


/**
 *  The following tests check if the default implementations for the RHF SCF solvers give the correct result.
 */


/**
 *  Check if our plain RHF SCF solver finds the correct energy. We will follow section 3.5.2 in Szabo.
 */
BOOST_AUTO_TEST_CASE(h2_sto3g_szabo_plain) {

    const double ref_total_energy = -1.1167;

    // Create the molecular Hamiltonian in an AO basis.
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {h2, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, h2);  // in an AO basis

    // Create a plain RHF SCF solver and solve the SCF equations.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    plain_rhf_scf_solver.perform(rhf_environment);


    // Check the total energy
    const double total_energy = rhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(h2).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-04);
}


/**
 *  Check if our plain RHF SCF solver finds results (energy, orbital energies and coefficient matrix) that are equal to results from HORTON.
 */
BOOST_AUTO_TEST_CASE(h2o_sto3g_horton_plain) {

    // List the reference data
    const double ref_total_energy = -74.942080055631;

    GQCP::VectorX<double> ref_orbital_energies {7};  // the STO-3G basisset has 7 basis functions for water
    ref_orbital_energies << -20.26289322, -1.20969863, -0.54796582, -0.43652631, -0.38758791, 0.47762043, 0.5881361;

    GQCP::SquareMatrix<double> ref_C_matrix {7};
    // clang-format off
    ref_C_matrix << -9.94434594e-01, -2.39158997e-01,  3.61117086e-17, -9.36837259e-02,  3.73303682e-31, -1.11639152e-01, -9.04958229e-17,
                    -2.40970260e-02,  8.85736467e-01, -1.62817254e-16,  4.79589270e-01, -1.93821120e-30,  6.69575233e-01,  5.16088339e-16,
                     1.59542752e-18,  5.29309704e-17, -6.07288675e-01, -1.49717339e-16,  8.94470461e-17, -8.85143477e-16,  9.19231270e-01,
                    -3.16155527e-03,  8.58957413e-02,  2.89059171e-16, -7.47426286e-01,  2.81871324e-30,  7.38494291e-01,  6.90314422e-16,
                     6.65079968e-35,  1.16150362e-32, -2.22044605e-16, -4.06685146e-30, -1.00000000e+00, -1.78495825e-31,  2.22044605e-16,
                     4.59373756e-03,  1.44038811e-01, -4.52995183e-01, -3.29475784e-01,  2.16823939e-16, -7.09847234e-01, -7.32462496e-01,
                     4.59373756e-03,  1.44038811e-01,  4.52995183e-01, -3.29475784e-01, -2.16823939e-16, -7.09847234e-01,  7.32462496e-01;
    // clang-format on
    const GQCP::RTransformation<double> ref_C {ref_C_matrix};

    // Perform our own RHF calculation.
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {water, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, water);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(water.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    plain_rhf_scf_solver.perform(rhf_environment);


    // Check the calculated results with the reference
    const double total_energy = rhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(water).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
    BOOST_CHECK(ref_orbital_energies.areEqualEigenvaluesAs(rhf_environment.orbital_energies.back(), 1.0e-06));
    BOOST_CHECK(ref_C.matrix().hasEqualSetsOfEigenvectorsAs(rhf_environment.coefficient_matrices.back().matrix(), 1.0e-05));
}


/**
 *  Check if the total RHF energy (calculated by our plain RHF SCF solver) for H2O matches the example from Crawdad. This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.
 */
BOOST_AUTO_TEST_CASE(crawdad_h2o_sto3g_plain) {

    const double ref_total_energy = -74.9420799281920;


    // Do our own RHF calculation
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
    // Check if the internuclear distance between O and H is really 1.1 A (= 2.07869 bohr), as specified in the text.
    BOOST_REQUIRE(std::abs(water.calculateInternuclearDistanceBetween(0, 1) - 2.07869) < 1.0e-4);

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {water, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, water);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(water.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    plain_rhf_scf_solver.perform(rhf_environment);


    // Check the total energy
    const double total_energy = rhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(water).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
}


/**
 *  Check if the total RHF energy for CH4 (calculated by our plain RHF SCF solver) matches the example from Crawdad. This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.
 */
BOOST_AUTO_TEST_CASE(crawdad_ch4_sto3g_plain) {

    const double ref_total_energy = -39.726850324347;

    // Do our own RHF calculation
    const auto methane = GQCP::Molecule::ReadXYZ("data/ch4_crawdad.xyz");
    // Check if the internuclear distance between C and H is really around 2.05 bohr, which is the bond distance Wikipedia (108.7 pm) specifies.
    BOOST_CHECK(std::abs(methane.calculateInternuclearDistanceBetween(0, 1) - 2.05) < 1.0e-1);

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {methane, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, methane);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(methane.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};


    plain_rhf_scf_solver.perform(rhf_environment);


    // Check the total energy
    const double total_energy = rhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(methane).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
}


/**
 *  Check if the total RHF energy for H2 (calculated by our plain RHF SCF solver) matches reference data from olsens (an implementation from Ayerslab). The reference data is for H2@RHF//6-31G** orbitals.
 */
BOOST_AUTO_TEST_CASE(h2_631gdp_plain) {

    const double ref_electronic_energy = -1.84444667247;

    // Perform our own RHF calculation.
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2_olsens.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {h2, "6-31G**"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, h2);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    plain_rhf_scf_solver.perform(rhf_environment);

    // Check the electronic energy
    BOOST_CHECK(std::abs(rhf_environment.electronic_energies.back() - ref_electronic_energy) < 1.0e-06);
}


/**
 *  Check if our damped RHF SCF solver finds results (energy, orbital energies and coefficient matrix) that are equal to results from HORTON.
 */
BOOST_AUTO_TEST_CASE(h2o_sto3g_horton_damped) {

    // List the reference data
    const double ref_total_energy = -74.942080055631;

    GQCP::VectorX<double> ref_orbital_energies {7};  // the STO-3G basisset has 7 basis functions for water
    ref_orbital_energies << -20.26289322, -1.20969863, -0.54796582, -0.43652631, -0.38758791, 0.47762043, 0.5881361;

    GQCP::SquareMatrix<double> ref_C_matrix {7};
    // clang-format off
    ref_C_matrix << -9.94434594e-01, -2.39158997e-01,  3.61117086e-17, -9.36837259e-02,  3.73303682e-31, -1.11639152e-01, -9.04958229e-17,
                    -2.40970260e-02,  8.85736467e-01, -1.62817254e-16,  4.79589270e-01, -1.93821120e-30,  6.69575233e-01,  5.16088339e-16,
                     1.59542752e-18,  5.29309704e-17, -6.07288675e-01, -1.49717339e-16,  8.94470461e-17, -8.85143477e-16,  9.19231270e-01,
                    -3.16155527e-03,  8.58957413e-02,  2.89059171e-16, -7.47426286e-01,  2.81871324e-30,  7.38494291e-01,  6.90314422e-16,
                     6.65079968e-35,  1.16150362e-32, -2.22044605e-16, -4.06685146e-30, -1.00000000e+00, -1.78495825e-31,  2.22044605e-16,
                     4.59373756e-03,  1.44038811e-01, -4.52995183e-01, -3.29475784e-01,  2.16823939e-16, -7.09847234e-01, -7.32462496e-01,
                     4.59373756e-03,  1.44038811e-01,  4.52995183e-01, -3.29475784e-01, -2.16823939e-16, -7.09847234e-01,  7.32462496e-01;
    // clang-format on
    const GQCP::RTransformation<double> ref_C {ref_C_matrix};

    // Perform our own RHF calculation.
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {water, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, water);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(water.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto damped_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DensityDamped(0.95);
    damped_rhf_scf_solver.perform(rhf_environment);


    // Check the calculated results with the reference
    const double total_energy = rhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(water).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
    BOOST_CHECK(ref_orbital_energies.areEqualEigenvaluesAs(rhf_environment.orbital_energies.back(), 1.0e-06));
    BOOST_CHECK(ref_C.matrix().hasEqualSetsOfEigenvectorsAs(rhf_environment.coefficient_matrices.back().matrix(), 1.0e-05));
}


/**
 *  Check if our DIIS RHF SCF solver finds the correct energy. We will follow section 3.5.2 in Szabo.
 */
BOOST_AUTO_TEST_CASE(h2_sto3g_szabo_diis) {

    const double ref_total_energy = -1.1167;

    // Create the molecular Hamiltonian in an AO basis.
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {h2, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, h2);  // in an AO basis

    // Create a DIIS RHF SCF solver and solve the SCF equations.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    diis_rhf_scf_solver.perform(rhf_environment);


    // Check the total energy
    const double total_energy = rhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(h2).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-04);
}


/**
 *  Check if our DIIS RHF SCF solver finds results (energy, orbital energies and coefficient matrix) that are equal to results from HORTON.
 */
BOOST_AUTO_TEST_CASE(h2o_sto3g_horton_diis) {

    // List the reference data
    const double ref_total_energy = -74.942080055631;

    GQCP::VectorX<double> ref_orbital_energies {7};  // the STO-3G basisset has 7 basis functions for water
    ref_orbital_energies << -20.26289322, -1.20969863, -0.54796582, -0.43652631, -0.38758791, 0.47762043, 0.5881361;

    GQCP::SquareMatrix<double> ref_C_matrix {7};
    // clang-format off
    ref_C_matrix << -9.94434594e-01, -2.39158997e-01,  3.61117086e-17, -9.36837259e-02,  3.73303682e-31, -1.11639152e-01, -9.04958229e-17,
                    -2.40970260e-02,  8.85736467e-01, -1.62817254e-16,  4.79589270e-01, -1.93821120e-30,  6.69575233e-01,  5.16088339e-16,
                     1.59542752e-18,  5.29309704e-17, -6.07288675e-01, -1.49717339e-16,  8.94470461e-17, -8.85143477e-16,  9.19231270e-01,
                    -3.16155527e-03,  8.58957413e-02,  2.89059171e-16, -7.47426286e-01,  2.81871324e-30,  7.38494291e-01,  6.90314422e-16,
                     6.65079968e-35,  1.16150362e-32, -2.22044605e-16, -4.06685146e-30, -1.00000000e+00, -1.78495825e-31,  2.22044605e-16,
                     4.59373756e-03,  1.44038811e-01, -4.52995183e-01, -3.29475784e-01,  2.16823939e-16, -7.09847234e-01, -7.32462496e-01,
                     4.59373756e-03,  1.44038811e-01,  4.52995183e-01, -3.29475784e-01, -2.16823939e-16, -7.09847234e-01,  7.32462496e-01;
    // clang-format on
    const GQCP::RTransformation<double> ref_C {ref_C_matrix};

    // perform our own RHF calculation.
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {water, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, water);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(water.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    diis_rhf_scf_solver.perform(rhf_environment);


    // Check the calculated results with the reference
    const double total_energy = rhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(water).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
    BOOST_CHECK(ref_orbital_energies.areEqualEigenvaluesAs(rhf_environment.orbital_energies.back(), 1.0e-06));
    BOOST_CHECK(ref_C.matrix().hasEqualSetsOfEigenvectorsAs(rhf_environment.coefficient_matrices.back().matrix(), 1.0e-05));
}


/**
 *  Check if the total RHF energy (calculated by our plain RHF SCF solver) for H2O matches the example from Crawdad. This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.
 */
BOOST_AUTO_TEST_CASE(crawdad_h2o_sto3g_diis) {

    const double ref_total_energy = -74.9420799281920;


    // Do our own RHF calculation
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
    // Check if the internuclear distance between O and H is really 1.1 A (= 2.07869 bohr), as specified in the text.
    BOOST_REQUIRE(std::abs(water.calculateInternuclearDistanceBetween(0, 1) - 2.07869) < 1.0e-4);

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {water, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, water);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(water.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    diis_rhf_scf_solver.perform(rhf_environment);


    // Check the total energy
    const double total_energy = rhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(water).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
}


/**
 *  Check if the total RHF energy for CH4 (calculated by our plain RHF SCF solver) matches the example from Crawdad. This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.
 */
BOOST_AUTO_TEST_CASE(crawdad_ch4_sto3g_diis) {

    const double ref_total_energy = -39.726850324347;


    // Do our own RHF calculation
    const auto methane = GQCP::Molecule::ReadXYZ("data/ch4_crawdad.xyz");
    // Check if the internuclear distance between C and H is really around 2.05 bohr, which is the bond distance Wikipedia (108.7 pm) specifies.
    BOOST_CHECK(std::abs(methane.calculateInternuclearDistanceBetween(0, 1) - 2.05) < 1.0e-1);

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {methane, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, methane);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(methane.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    diis_rhf_scf_solver.perform(rhf_environment);


    // Check the total energy
    const double total_energy = rhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(methane).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
}


/**
 *  Check if the total RHF energy for H2 (calculated by our plain RHF SCF solver) matches reference data from olsens (an implementation from Ayerslab). The reference data is for H2@RHF//6-31G** orbitals.
 */
BOOST_AUTO_TEST_CASE(h2_631gdp_diis) {

    const double ref_electronic_energy = -1.84444667247;

    // Perform our own RHF calculation.
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2_olsens.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {h2, "6-31G**"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, h2);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    diis_rhf_scf_solver.perform(rhf_environment);


    // Check the electronic energy
    BOOST_CHECK(std::abs(rhf_environment.electronic_energies.back() - ref_electronic_energy) < 1.0e-06);
}
