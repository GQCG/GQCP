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

#define BOOST_TEST_MODULE "SpinResolvedONV"

#include <boost/test/unit_test.hpp>

#include "ONVBasis/SpinResolvedONV.hpp"


/**
 *  Check if the creation of the 'RHF' ONV works as expected.
 */
BOOST_AUTO_TEST_CASE(RHF) {

    // For K=5 alpha-spinorbitals, K=5 beta-spinorbitals, and N_P=3 electron electron pairs, the 'RHF' ONV should be ("00111" = 7) ("00111" = 7)
    const size_t K = 5;
    const size_t N_P = 3;
    const GQCP::SpinUnresolvedONV alpha_onv {K, N_P, 7};

    const GQCP::SpinResolvedONV reference {alpha_onv, alpha_onv};
    BOOST_CHECK(reference == GQCP::SpinResolvedONV::RHF(K, N_P));
}


/**
 *  Check if the creation of the 'UHF' ONV works as expected.
 */
BOOST_AUTO_TEST_CASE(UHF) {

    // For K=5 alpha-spinorbitals, K=5 beta-spinorbitals, and N_alpha=3 alpha-electrons and N_beta=2 beta-electrons, the 'UHF' ONV should be ("00111" = 7) ("00011" = 3)
    const size_t K = 5;
    const size_t N_alpha = 3;
    const size_t N_beta = 2;
    const GQCP::SpinUnresolvedONV alpha_onv {K, N_alpha, 7};
    const GQCP::SpinUnresolvedONV beta_onv {K, N_beta, 3};

    const GQCP::SpinResolvedONV reference {alpha_onv, beta_onv};
    BOOST_CHECK(reference == GQCP::SpinResolvedONV::UHF(K, N_alpha, N_beta));
}
