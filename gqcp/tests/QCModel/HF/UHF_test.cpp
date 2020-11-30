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

#define BOOST_TEST_MODULE "QCModel::UHF"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/UHF/UHF.hpp"
#include "QCMethod/HF/UHF/UHFSCFSolver.hpp"
#include "QCModel/HF/UHF.hpp"


/**
 *  Check if the UHF basic constructor throws as expected.
 */
BOOST_AUTO_TEST_CASE(basic_constructor) {

    // Numbers of electrons.
    const size_t N_alpha1 = 1;
    const size_t N_alpha2 = 2;
    const size_t N_beta1 = 1;
    const size_t N_beta2 = 2;

    // Numbers of spatial orbitals.
    const size_t K_alpha1 = 5;
    const size_t K_alpha2 = 1;
    const size_t K_beta1 = 5;
    const size_t K_beta2 = 1;

    // Orbital energies.
    const GQCP::VectorX<double> orbital_energies_alpha1 = GQCP::VectorX<double>::Zero(K_alpha1);
    const GQCP::VectorX<double> orbital_energies_beta1 = GQCP::VectorX<double>::Zero(K_beta1);
    const GQCP::VectorX<double> orbital_energies_alpha2 = GQCP::VectorX<double>::Zero(K_alpha2);
    const GQCP::VectorX<double> orbital_energies_beta2 = GQCP::VectorX<double>::Zero(K_beta2);

    // Transformation matrices.
    const GQCP::UTransformationComponent<double> C_alpha1 = GQCP::UTransformationComponent<double>::Identity(K_alpha1);
    const GQCP::UTransformationComponent<double> C_beta1 = GQCP::UTransformationComponent<double>::Identity(K_beta1);
    const GQCP::UTransformationComponent<double> C_alpha2 = GQCP::UTransformationComponent<double>::Identity(K_alpha2);
    const GQCP::UTransformationComponent<double> C_beta2 = GQCP::UTransformationComponent<double>::Identity(K_beta2);

    const GQCP::UTransformation<double> C_1 = GQCP::UTransformation<double> {C_alpha1, C_beta1};
    const GQCP::UTransformation<double> C_2 = GQCP::UTransformation<double> {C_alpha2, C_beta2};
    const GQCP::UTransformation<double> C_3 = GQCP::UTransformation<double> {C_alpha1, C_beta2};
    const GQCP::UTransformation<double> C_4 = GQCP::UTransformation<double> {C_alpha2, C_beta1};


    // Check if an expected correct constructor doesn't throw.
    BOOST_CHECK_NO_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta1, orbital_energies_alpha1, orbital_energies_beta1, C_1));

    // Check throws if there are too many electrons.
    BOOST_CHECK_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta2, orbital_energies_alpha2, orbital_energies_beta2, C_2), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::QCModel::UHF<double>(N_alpha2, N_beta1, orbital_energies_alpha2, orbital_energies_beta2, C_2), std::invalid_argument);

    // Check throws if the number of spatial orbitals (dimension of C) doesn't match the number of orbital energies.
    BOOST_CHECK_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta1, orbital_energies_alpha1, orbital_energies_beta1, C_4), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta1, orbital_energies_alpha1, orbital_energies_beta1, C_3), std::invalid_argument);


    // Check that the constructor that sets the orbital energies to zeros doesn't throw.
    BOOST_CHECK_NO_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta1, C_1));
}


/**
 *  Check if the methods for returning spin-orbital energies are correctly implemented.
 */
BOOST_AUTO_TEST_CASE(spinorbitalEnergies) {

    // Set up toy UHF model parameters.
    const size_t K = 2;
    const auto C = GQCP::RTransformation<double>::Identity(K);
    GQCP::VectorX<double> orbital_energies {K};
    orbital_energies << -0.5, 0.5;

    GQCP::QCModel::RHF<double> rhf_parameters {1, orbital_energies, C};
    GQCP::QCModel::UHF<double> uhf_parameters {rhf_parameters};


    // Provide reference values and check the results.
    GQCP::VectorX<double> ref_spinorbital_energies_blocked {2 * K};
    ref_spinorbital_energies_blocked << -0.5, 0.5, -0.5, 0.5;
    BOOST_CHECK(uhf_parameters.spinOrbitalEnergiesBlocked().isApprox(ref_spinorbital_energies_blocked, 1.0e-12));
}


/**
 *  Check if the UHF energy is equal to the expectation value of the Hamiltonian through its density matrices.
 */
BOOST_AUTO_TEST_CASE(UHF_DMs) {

    // Perform a UHF calculation.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto r_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In an AO basis.

    const auto N_a = molecule.numberOfElectronPairs();  // The number of alpha electrons.
    const auto N_b = molecule.numberOfElectronPairs();  // The number of beta electrons.
    auto uhf_environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_a, N_b, r_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_uhf_scf_solver = GQCP::UHFSCFSolver<double>::Plain();

    const auto uhf_qc_structure = GQCP::QCMethod::UHF<double>().optimize(plain_uhf_scf_solver, uhf_environment);
    const auto uhf_parameters = uhf_qc_structure.groundStateParameters();
    const auto uhf_energy = uhf_qc_structure.groundStateEnergy();

    // Determine the RHF energy through the expectation value of the Hamiltonian, and check the result.
    // Do the calculations in the RHF MO basis, in order to check the implementation of the RHF density matrices in MO basis.
    auto u_spin_orbital_basis = GQCP::USpinOrbitalBasis<double, GQCP::GTOShell>::FromRestricted(spin_orbital_basis);
    u_spin_orbital_basis.transform(uhf_parameters.expansion());
    auto u_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(u_spin_orbital_basis, molecule);

    const auto D_MO = uhf_parameters.calculateOrthonormalBasis1DM();
    const auto d_MO = uhf_parameters.calculateOrthonormalBasis2DM();
    const double expectation_value = u_hamiltonian.calculateExpectationValue(D_MO, d_MO);

    BOOST_CHECK(std::abs(uhf_energy - expectation_value) < 1.0e-12);
}
