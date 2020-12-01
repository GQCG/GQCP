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
BOOST_AUTO_TEST_CASE(dyson_amplitudes) {

    // Set up the manually calculated references.
    GQCP::Vector<double, 2> reference_amplitudes_alpha {0.39739531532399996, 0.9116729926689999};
    GQCP::Vector<double, 2> reference_amplitudes_beta = {0.537653264399, 0.794791398869};


    // Set up the toy linear expansions.
    const size_t K = 2;
    const size_t N = 2;

    GQCP::VectorX<double> coeffs1 = GQCP::VectorX<double>::Zero(4);
    coeffs1 << 0.182574, 0.365148, 0.547723, 0.730297;
    GQCP::VectorX<double> coeffs2 = GQCP::VectorX<double>::Zero(2);
    coeffs2 << 0.640184, 0.768221;

    const GQCP::SpinResolvedONVBasis onv_basis {K, N / 2, N / 2};            // The reference ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis_alpha {K, N / 2 - 1, N / 2};  // An ONV basis with one less alpha electron.
    const GQCP::SpinResolvedONVBasis onv_basis_beta {K, N / 2, N / 2 - 1};   // An ONV basis with one less beta electron.

    const auto linear_expansion = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>(onv_basis, coeffs1);
    const auto linear_expansion_alpha = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>(onv_basis_alpha, coeffs2);
    const auto linear_expansion_beta = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>(onv_basis_beta, coeffs2);


    // Calculate the Dyson amplitudes for both situations (alpha-reference) and (beta-reference), and check with the manual calculations.
    const auto dyson_orbital_alpha = GQCP::DysonOrbital<double>::TransitionAmplitudes(linear_expansion, linear_expansion_alpha);
    const auto& dyson_coefficients_alpha = dyson_orbital_alpha.amplitudes();

    const auto dyson_orbital_beta = GQCP::DysonOrbital<double>::TransitionAmplitudes(linear_expansion, linear_expansion_beta);
    const auto& dyson_coefficients_beta = dyson_orbital_beta.amplitudes();

    BOOST_CHECK(dyson_coefficients_alpha.isApprox(reference_amplitudes_alpha, 1.0e-06));
    BOOST_CHECK(dyson_coefficients_beta.isApprox(reference_amplitudes_beta, 1.0e-06));
}
