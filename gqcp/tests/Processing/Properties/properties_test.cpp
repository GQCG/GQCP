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
#include "Operator/FirstQuantized/Operator.hpp"
#include "Processing/Properties/RHFElectricalResponseSolver.hpp"
#include "Processing/Properties/properties.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "Utilities/units.hpp"


/**
 *  Check the calculation of the CO dipole moment from a CCCBDB reference value.
 */
BOOST_AUTO_TEST_CASE(dipole_CO_STO_3G) {

    // Initialize the molecule and molecular Hamiltonian for CO.
    const GQCP::Nucleus C {6, 0.0, 0.0, 0.0};
    const GQCP::Nucleus O {8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145)};  // From CCCBDB, STO-3G geometry.
    const GQCP::Molecule molecule {{C, O}};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.

    const auto K = spin_orbital_basis.numberOfSpatialOrbitals();
    const size_t N = molecule.numberOfElectrons();

    // Solve the RHF SCF equations.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    const double total_energy = rhf_qc_structure.groundStateEnergy() + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_REQUIRE(std::abs(total_energy - (-111.225)) < 1.0e-02);  // From CCCBDB, require a correct RHF solution to be found.


    // Calculate the RHF 1-DM and the dipole operator in RHF MO basis.
    const auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);
    auto dipole_op = spin_orbital_basis.quantize(GQCP::Operator::ElectronicDipole());
    dipole_op.transform(rhf_parameters.expansion());

    // Calculate the RHF total dipole moment in the MO basis and check with the reference value.
    GQCP::Vector<double, 3> total_dipole_moment = GQCP::Operator::NuclearDipole(molecule).value() + dipole_op.calculateExpectationValue(D).asVector();
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
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.

    const auto K = spin_orbital_basis.numberOfSpatialOrbitals();
    const auto N = molecule.numberOfElectrons();

    // Solve the RHF SCF equations.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    const double total_energy = rhf_qc_structure.groundStateEnergy() + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_REQUIRE(std::abs(total_energy - (-107.500654)) < 1.0e-05);  // From CCCBDB, require a correct RHF solution to be found.


    // Calculate the RHF 1-DM and the dipole operator in RHF MO basis.
    const auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);
    auto dipole_op = spin_orbital_basis.quantize(GQCP::Operator::ElectronicDipole());
    dipole_op.transform(rhf_parameters.expansion());

    // Calculate the RHF total dipole moment in the MO basis and check with the reference value.
    GQCP::Vector<double, 3> total_dipole_moment = GQCP::Operator::NuclearDipole(molecule).value() + dipole_op.calculateExpectationValue(D).asVector();
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
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.


    // Do the RHF calculation to get the canonical RHF orbitals.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective(sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();


    // Transform the orbitals to the RHF basis and prepare the dipole integrals in the RHF basis.
    GQCP::transform(rhf_parameters.expansion(), spin_orbital_basis, sq_hamiltonian);
    const auto dipole_op = spin_orbital_basis.quantize(GQCP::Operator::ElectronicDipole());


    // Find the RHF wave function response.
    GQCP::RHFElectricalResponseSolver cphf_solver {molecule.numberOfElectrons() / 2};
    const auto x = cphf_solver.calculateWaveFunctionResponse(sq_hamiltonian, dipole_op);


    // Calculate the RHF polarizability and check with the reference value.
    const auto F_p = cphf_solver.calculateParameterResponseForce(dipole_op);
    const auto alpha = GQCP::calculateElectricPolarizability(F_p, x);
    const auto alpha_zz = alpha(2, 2);

    BOOST_CHECK(std::abs(alpha_zz - ref_alpha_zz) < 1.0e-05);
}
