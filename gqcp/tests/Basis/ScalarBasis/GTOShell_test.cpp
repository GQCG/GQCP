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

#define BOOST_TEST_MODULE "GTOShell"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis/GTOShell.hpp"


/**
 *  Check if the GTOShell constructor throws as expected.
 */
BOOST_AUTO_TEST_CASE(Shell_constructor_throws) {


    const std::vector<double> exp1 {1.0, 1.1};
    const std::vector<double> coeff1 {0.5, 1.0};
    const std::vector<double> coeff2 {0.5, 1.0, 1.5};

    const GQCP::GTOShell shell1 {0, GQCP::Nucleus(), exp1, coeff1};                                     // shouldn't throw
    BOOST_CHECK_THROW(GQCP::GTOShell shell2(0, GQCP::Nucleus(), exp1, coeff2), std::invalid_argument);  // too many coefficients for the number of exponents
}


/**
 *  Check if the number of basis functions is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(numberOfBasisFunctions) {

    const std::vector<double> exp {1.0, 1.1};
    const std::vector<double> coeff {0.5, 1.0};
    const GQCP::Nucleus nucleus {};


    // Check the Cartesian shell implementation.
    const GQCP::GTOShell cartesian_s_shell {0, nucleus, exp, coeff, false};
    const GQCP::GTOShell cartesian_p_shell {1, nucleus, exp, coeff, false};
    const GQCP::GTOShell cartesian_d_shell {2, nucleus, exp, coeff, false};
    const GQCP::GTOShell cartesian_f_shell {3, nucleus, exp, coeff, false};

    BOOST_CHECK(cartesian_s_shell.numberOfBasisFunctions() == 1);
    BOOST_CHECK(cartesian_p_shell.numberOfBasisFunctions() == 3);
    BOOST_CHECK(cartesian_d_shell.numberOfBasisFunctions() == 6);
    BOOST_CHECK(cartesian_f_shell.numberOfBasisFunctions() == 10);


    // Check the spherical shell implementation.
    const GQCP::GTOShell spherical_s_shell {0, nucleus, exp, coeff};
    const GQCP::GTOShell spherical_p_shell {1, nucleus, exp, coeff};
    const GQCP::GTOShell spherical_d_shell {2, nucleus, exp, coeff};
    const GQCP::GTOShell spherical_f_shell {3, nucleus, exp, coeff};

    BOOST_CHECK(spherical_s_shell.numberOfBasisFunctions() == 1);
    BOOST_CHECK(spherical_p_shell.numberOfBasisFunctions() == 3);
    BOOST_CHECK(spherical_d_shell.numberOfBasisFunctions() == 5);
    BOOST_CHECK(spherical_f_shell.numberOfBasisFunctions() == 7);
}


/**
 *  Check if operator== is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(operator_equals) {

    const size_t l1 = 0;
    const size_t l2 = 1;
    const GQCP::Nucleus nucleus1 {};
    const GQCP::Nucleus nucleus2 {1, 0.0, 0.0, 0.0};
    const std::vector<double> exp1 {1.0, 1.1};
    const std::vector<double> exp2 {2.0, 2.1};
    const std::vector<double> coeff1 {0.5, 1.0};
    const std::vector<double> coeff2 {2.5, 3.0};


    // Test for equality.
    BOOST_CHECK(GQCP::GTOShell(l1, nucleus1, exp1, coeff1) == GQCP::GTOShell(l1, nucleus1, exp1, coeff1));

    // Check if a different angular momentum causes inequality.
    BOOST_CHECK(!(GQCP::GTOShell(l1, nucleus1, exp1, coeff1) == GQCP::GTOShell(l2, nucleus1, exp1, coeff1)));

    // Check if a different nucleus causes inequality.
    BOOST_CHECK(!(GQCP::GTOShell(l1, nucleus1, exp1, coeff1) == GQCP::GTOShell(l1, nucleus2, exp1, coeff1)));

    // Check if different Gaussian exponents cause inequality.
    BOOST_CHECK(!(GQCP::GTOShell(l1, nucleus1, exp1, coeff1) == GQCP::GTOShell(l1, nucleus1, exp2, coeff1)));

    // Check if different contraction coefficients cause inequality.
    BOOST_CHECK(!(GQCP::GTOShell(l1, nucleus1, exp1, coeff1) == GQCP::GTOShell(l1, nucleus1, exp1, coeff2)));
}


/**
 *  Check if embedding and unembedding the normalization factors of the individual primitives works as expected.
 */
BOOST_AUTO_TEST_CASE(embed_normalization_factor_primitives) {

    const std::vector<double> exp {1.0, 1.1};
    const GQCP::Nucleus nucleus {};
    GQCP::GTOShell s_shell {0, nucleus, exp, {0.5, 1.0}};
    auto s_shell_copy = s_shell;

    const GQCP::GTOShell ref_embedded_shell {0, nucleus, exp, {0.3563527351774951, 0.7655165883890895}};  // manual calculation

    // Check if embedding and un-embedding is a zero operation.
    // Also test the behavior of calling the normalization functions twice.
    s_shell_copy.embedNormalizationFactorsOfPrimitives();
    s_shell_copy.embedNormalizationFactorsOfPrimitives();
    BOOST_CHECK(s_shell_copy == ref_embedded_shell);

    s_shell_copy.unEmbedNormalizationFactorsOfPrimitives();
    BOOST_CHECK(s_shell == s_shell_copy);

    s_shell_copy.unEmbedNormalizationFactorsOfPrimitives();
    BOOST_CHECK(s_shell == s_shell_copy);
}


/**
 *  Check if embedding the total normalization factor of the contraction works as expected.
 */
BOOST_AUTO_TEST_CASE(embed_total_normalization_factor) {

    const std::vector<double> exp {1.0, 1.1};
    const GQCP::Nucleus nucleus {};
    GQCP::GTOShell s_shell {0, nucleus, exp, {0.5, 1.0}};

    const GQCP::GTOShell ref_embedded_shell {0, nucleus, exp, {0.3854188481329033, 0.7708376962658066}};  // manual calculation


    // Check if the embedding of the total normalization factor is correct.
    s_shell.embedNormalizationFactor();
    BOOST_CHECK(ref_embedded_shell == s_shell);
}
