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

#define BOOST_TEST_MODULE "DysonOrbital"

#include <boost/test/unit_test.hpp>

#include "Processing/Properties/DysonOrbital.hpp"


/**
 *  Test the algorithm for the Dyson amplitudes against manually calculated coefficients for two (normalized) toy wave functions.
 */
BOOST_AUTO_TEST_CASE(dyson_amplitudes_spin_resolved_1) {

    // Set up the manually calculated references.
    GQCP::Vector<double, 2> reference_amplitudes_alpha {0.39739531532399996, 0.9116729926689999};
    GQCP::Vector<double, 2> reference_amplitudes_beta = {0.537653264399, 0.794791398869};


    // Set up the toy linear expansions.
    const size_t K = 2;
    const size_t N = 2;

    GQCP::VectorX<double> coeffs_J = GQCP::VectorX<double>::Zero(4);
    coeffs_J << 0.182574, 0.365148, 0.547723, 0.730297;
    GQCP::VectorX<double> coeffs_I = GQCP::VectorX<double>::Zero(2);
    coeffs_I << 0.640184, 0.768221;

    const GQCP::SpinResolvedONVBasis onv_basis {K, N / 2, N / 2};            // The reference ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis_alpha {K, N / 2 - 1, N / 2};  // An ONV basis with one less alpha electron.
    const GQCP::SpinResolvedONVBasis onv_basis_beta {K, N / 2, N / 2 - 1};   // An ONV basis with one less beta electron.

    const auto linear_expansion = GQCP::LinearExpansion<double, GQCP::SpinResolvedONVBasis>(onv_basis, coeffs_J);
    const auto linear_expansion_alpha = GQCP::LinearExpansion<double, GQCP::SpinResolvedONVBasis>(onv_basis_alpha, coeffs_I);
    const auto linear_expansion_beta = GQCP::LinearExpansion<double, GQCP::SpinResolvedONVBasis>(onv_basis_beta, coeffs_I);


    // Calculate the Dyson amplitudes for both situations (alpha-reference) and (beta-reference), and check with the manual calculations.
    const auto dyson_orbital_alpha = GQCP::DysonOrbital<double>::TransitionAmplitudes(linear_expansion, linear_expansion_alpha);
    const auto& dyson_coefficients_alpha = dyson_orbital_alpha.amplitudes();

    const auto dyson_orbital_beta = GQCP::DysonOrbital<double>::TransitionAmplitudes(linear_expansion, linear_expansion_beta);
    const auto& dyson_coefficients_beta = dyson_orbital_beta.amplitudes();

    BOOST_CHECK(dyson_coefficients_alpha.isApprox(reference_amplitudes_alpha, 1.0e-06));
    BOOST_CHECK(dyson_coefficients_beta.isApprox(reference_amplitudes_beta, 1.0e-06));
}

/**
 *  Test the algorithm for the Dyson amplitudes against manually calculated coefficients for two (normalized) toy wave functions where each spin component has more than one electron. This is to check whether the sign of the Dyson amplitudes are correct.
 */
BOOST_AUTO_TEST_CASE(dyson_amplitudes_spin_resolved_2) {

    // Set up the manually calculated references.
    GQCP::Vector<double, 3> reference_amplitudes_alpha {0.578739438503937, -0.11202006721497709, -0.8605609287217918};
    GQCP::Vector<double, 3> reference_amplitudes_beta = {0.8606702596667963, 0.0868282794535802, -0.8722314733220338};


    // Set up the toy linear expansions.
    const size_t K = 3;
    const size_t N = 4;

    GQCP::VectorX<double> coeffs_J = GQCP::VectorX<double>::Zero(9);
    coeffs_J << 0.56494513, 0.38187498, 0.82585997, 0.23923204, 0.32256349, 0.22982795, 0.22972143, 0.71964626, 0.41650422;
    GQCP::VectorX<double> coeffs_I = GQCP::VectorX<double>::Zero(9);
    coeffs_I << 0.26297864, 0.47549281, 0.13096657, 0.11171302, 0.79625911, 0.03717573, 0.12209374, 0.17338261, 0.4164798;

    const GQCP::SpinResolvedONVBasis onv_basis {K, N / 2, N / 2};            // The reference ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis_alpha {K, N / 2 - 1, N / 2};  // An ONV basis with one less alpha electron.
    const GQCP::SpinResolvedONVBasis onv_basis_beta {K, N / 2, N / 2 - 1};   // An ONV basis with one less beta electron.

    const auto linear_expansion = GQCP::LinearExpansion<double, GQCP::SpinResolvedONVBasis>(onv_basis, coeffs_J);
    const auto linear_expansion_alpha = GQCP::LinearExpansion<double, GQCP::SpinResolvedONVBasis>(onv_basis_alpha, coeffs_I);
    const auto linear_expansion_beta = GQCP::LinearExpansion<double, GQCP::SpinResolvedONVBasis>(onv_basis_beta, coeffs_I);


    // Calculate the Dyson amplitudes for both situations (alpha-reference) and (beta-reference), and check with the manual calculations.
    const auto dyson_orbital_alpha = GQCP::DysonOrbital<double>::TransitionAmplitudes(linear_expansion, linear_expansion_alpha);
    const auto& dyson_coefficients_alpha = dyson_orbital_alpha.amplitudes();
    const auto dyson_orbital_beta = GQCP::DysonOrbital<double>::TransitionAmplitudes(linear_expansion, linear_expansion_beta);
    const auto& dyson_coefficients_beta = dyson_orbital_beta.amplitudes();

    BOOST_CHECK(dyson_coefficients_alpha.isApprox(reference_amplitudes_alpha, 1.0e-06));
    BOOST_CHECK(dyson_coefficients_beta.isApprox(reference_amplitudes_beta, 1.0e-06));
}

/**
 *  Test the algorithm for the Dyson amplitudes against manually calculated coefficients for two (normalized) toy wave functions with more than one electron.
 */
BOOST_AUTO_TEST_CASE(dyson_amplitudes_spin_unresolved) {

    // Set up the manually calculated references.
    GQCP::Vector<double, 3> reference_amplitudes {0.49666958935692024, 0.332773392216452, -0.7910269806599591};

    // Set up the toy linear expansions.
    const size_t K = 3;
    const size_t N = 2;

    GQCP::VectorX<double> coeffs_I = GQCP::VectorX<double>::Zero(3);
    coeffs_I << 0.89044211, 0.36096181, 0.71094162;
    GQCP::VectorX<double> coeffs_J = GQCP::VectorX<double>::Zero(3);
    coeffs_J << 0.11975192, 0.63780725, 0.61806136;

    const GQCP::SpinUnresolvedONVBasis onv_basis_J {K, N};      // The reference ONV basis.
    const GQCP::SpinUnresolvedONVBasis onv_basis_I {K, N - 1};  // An ONV basis with one less electron.

    const auto linear_expansion_J = GQCP::LinearExpansion<double, GQCP::SpinUnresolvedONVBasis>(onv_basis_J, coeffs_J);
    const auto linear_expansion_I = GQCP::LinearExpansion<double, GQCP::SpinUnresolvedONVBasis>(onv_basis_I, coeffs_I);

    // Calculate the Dyson orbital and check with manual calculations.
    const auto dyson_orbital = GQCP::DysonOrbital<double>::TransitionAmplitudes(linear_expansion_J, linear_expansion_I);
    const auto& dyson_coefficients = dyson_orbital.amplitudes();

    BOOST_CHECK(dyson_coefficients.isApprox(reference_amplitudes));
}
