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
    const GQCP::TransformationMatrix<double> C_alpha1 = GQCP::TransformationMatrix<double>::Identity(K_alpha1, K_alpha1);
    const GQCP::TransformationMatrix<double> C_beta1 = GQCP::TransformationMatrix<double>::Identity(K_beta1, K_beta1);
    const GQCP::TransformationMatrix<double> C_alpha2 = GQCP::TransformationMatrix<double>::Identity(K_alpha2, K_alpha2);
    const GQCP::TransformationMatrix<double> C_beta2 = GQCP::TransformationMatrix<double>::Identity(K_beta2, K_beta2);


    // Check if an expected correct constructor doesn't throw.
    BOOST_CHECK_NO_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta1, orbital_energies_alpha1, orbital_energies_beta1, C_alpha1, C_beta1));

    // Check throws if there are too many electrons.
    BOOST_CHECK_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta2, orbital_energies_alpha2, orbital_energies_beta2, C_alpha2, C_beta2), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::QCModel::UHF<double>(N_alpha2, N_beta1, orbital_energies_alpha2, orbital_energies_beta2, C_alpha2, C_beta2), std::invalid_argument);

    // Check throws if the number of spatial orbitals (dimension of C) doesn't match the number of orbital energies.
    BOOST_CHECK_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta1, orbital_energies_alpha1, orbital_energies_beta1, C_alpha2, C_beta1), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta1, orbital_energies_alpha1, orbital_energies_beta1, C_alpha1, C_beta2), std::invalid_argument);


    // Check that the constructor that sets the orbital energies to zeros doesn't throw.
    BOOST_CHECK_NO_THROW(GQCP::QCModel::UHF<double>(N_alpha1, N_beta1, C_alpha1, C_beta1));
}


/**
 *  Check if the methods for returning spinorbital energies are correctly implemented.
 */
BOOST_AUTO_TEST_CASE(spinorbitalEnergies) {

    // Set up toy UHF model parameters.
    const size_t K = 2;
    const GQCP::TransformationMatrix<double> C = GQCP::TransformationMatrix<double>::Identity(K, K);
    GQCP::VectorX<double> orbital_energies {K};
    orbital_energies << -0.5, 0.5;

    GQCP::QCModel::RHF<double> rhf_parameters {1, orbital_energies, C};
    GQCP::QCModel::UHF<double> uhf_parameters {rhf_parameters};


    // Provide reference values and check the results.
    GQCP::VectorX<double> ref_spinorbital_energies_blocked {2 * K};
    ref_spinorbital_energies_blocked << -0.5, 0.5, -0.5, 0.5;
    BOOST_CHECK(uhf_parameters.spinOrbitalEnergiesBlocked().isApprox(ref_spinorbital_energies_blocked, 1.0e-12));
}
