// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#define BOOST_TEST_MODULE "Shell"

#include <boost/test/unit_test.hpp>

#include "Basis/Shell.hpp"


BOOST_AUTO_TEST_CASE ( Shell_constructor_throws ) {

    std::vector<double> exp1 {1.0, 1.1};
    std::vector<double> coeff1 {0.5, 1.0};
    std::vector<double> coeff2 {0.5, 1.0, 1.5};

    GQCP::Shell shell1 (0, GQCP::Atom(), exp1, coeff1);
    BOOST_CHECK_THROW(GQCP::Shell shell2 (0, GQCP::Atom(), exp1, coeff2), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( numberOfBasisFunctions ) {

    std::vector<double> exp {1.0, 1.1};
    std::vector<double> coeff {0.5, 1.0};
    GQCP::Atom atom {};
    GQCP::Shell cartesian_s_shell (0, atom, exp, coeff, false);
    GQCP::Shell cartesian_p_shell (1, atom, exp, coeff, false);
    GQCP::Shell cartesian_d_shell (2, atom, exp, coeff, false);
    GQCP::Shell cartesian_f_shell (3, atom, exp, coeff, false);

    BOOST_CHECK(cartesian_s_shell.numberOfBasisFunctions() == 1);
    BOOST_CHECK(cartesian_p_shell.numberOfBasisFunctions() == 3);
    BOOST_CHECK(cartesian_d_shell.numberOfBasisFunctions() == 6);
    BOOST_CHECK(cartesian_f_shell.numberOfBasisFunctions() == 10);


    GQCP::Shell spherical_s_shell (0, atom, exp, coeff);
    GQCP::Shell spherical_p_shell (1, atom, exp, coeff);
    GQCP::Shell spherical_d_shell (2, atom, exp, coeff);
    GQCP::Shell spherical_f_shell (3, atom, exp, coeff);

    BOOST_CHECK(spherical_s_shell.numberOfBasisFunctions() == 1);
    BOOST_CHECK(spherical_p_shell.numberOfBasisFunctions() == 3);
    BOOST_CHECK(spherical_d_shell.numberOfBasisFunctions() == 5);
    BOOST_CHECK(spherical_f_shell.numberOfBasisFunctions() == 7);
}


BOOST_AUTO_TEST_CASE ( operator_equals ) {

    size_t l1 = 0;
    size_t l2 = 1;
    GQCP::Atom atom1 {};
    GQCP::Atom atom2 (1.0,  0.0, 0.0, 0.0);
    std::vector<double> exp1 {1.0, 1.1};
    std::vector<double> exp2 {2.0, 2.1};
    std::vector<double> coeff1 {0.5, 1.0};
    std::vector<double> coeff2 {2.5, 3.0};


    // Test for equality
    BOOST_CHECK(GQCP::Shell(l1, atom1, exp1, coeff1) == GQCP::Shell(l1, atom1, exp1, coeff1));

    // Check if a different angular momentum causes inequality
    BOOST_CHECK(!(GQCP::Shell(l1, atom1, exp1, coeff1) == GQCP::Shell(l2, atom1, exp1, coeff1)));

    // Check if a different atom causes inequality
    BOOST_CHECK(!(GQCP::Shell(l1, atom1, exp1, coeff1) == GQCP::Shell(l1, atom2, exp1, coeff1)));

    // Check if different Gaussian exponents cause inequality
    BOOST_CHECK(!(GQCP::Shell(l1, atom1, exp1, coeff1) == GQCP::Shell(l1, atom1, exp2, coeff1)));

    // Check if different contraction coefficients cause inequality
    BOOST_CHECK(!(GQCP::Shell(l1, atom1, exp1, coeff1) == GQCP::Shell(l1, atom1, exp1, coeff2)));
}


BOOST_AUTO_TEST_CASE ( embed_normalization_factor_primitives ) {

    std::vector<double> exp {1.0, 1.1};
    GQCP::Atom atom {};
    GQCP::Shell s_shell (0, atom, exp, {0.5, 1.0});
    auto s_shell_copy = s_shell;

    GQCP::Shell ref_embedded_shell (0, atom, exp, {0.3563527351774951, 0.7655165883890895});  // manual calculation

    // Check if embedding and un-embedding is a zero operation
    // Also test the behavior of calling the normalization functions twice
    s_shell_copy.embedNormalizationFactorsOfPrimitives();
    s_shell_copy.embedNormalizationFactorsOfPrimitives();
    BOOST_CHECK(s_shell_copy == ref_embedded_shell);

    s_shell_copy.unEmbedNormalizationFactorsOfPrimitives();
    BOOST_CHECK(s_shell == s_shell_copy);

    s_shell_copy.unEmbedNormalizationFactorsOfPrimitives();
    BOOST_CHECK(s_shell == s_shell_copy);
}


BOOST_AUTO_TEST_CASE ( embed_total_normalization_factor ) {

    std::vector<double> exp {1.0, 1.1};
    GQCP::Atom atom {};
    GQCP::Shell s_shell (0, atom, exp, {0.5, 1.0});

    GQCP::Shell ref_embedded_shell (0, atom, exp, {0.3854188481329033, 0.7708376962658066});  // manual calculation


    // Check if the embedding of the total normalization factor is correct
    s_shell.embedNormalizationFactor();
    BOOST_CHECK(ref_embedded_shell == s_shell);
}
