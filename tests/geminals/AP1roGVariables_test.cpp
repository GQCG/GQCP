// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2018  the GQCG developers
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
#define BOOST_TEST_MODULE "AP1roGVariables"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "geminals/AP1roGVariables.hpp"


BOOST_AUTO_TEST_CASE ( constructor ) {

    // Check a correct constructor
    GQCP::AP1roGVariables x (4, 6);

    // We can't create 4 geminals in 4 orbitals
    BOOST_CHECK_THROW(GQCP::AP1roGVariables (4, 4), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( constructor_vector ) {

    // Check a correct constructor
    Eigen::VectorXd x = Eigen::VectorXd::Zero(6);
    BOOST_CHECK_NO_THROW(GQCP::AP1roGVariables (x, 2, 5));

    // Check wrong parameters N_P and K
    BOOST_CHECK_THROW(GQCP::AP1roGVariables (x, 2, 4), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::AP1roGVariables (x, 1, 5), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( asMatrix ) {

    // For N_P=2 and K=5, we have an AP1roG variables matrix that looks like the following matrix:
    Eigen::MatrixXd X (2, 3);
    X << 1, 2, 3,
         4, 5, 6;

    // The variables arranged in a vector are then represented by:
    Eigen::VectorXd x (6);
    x << 1, 2, 3, 4, 5, 6;


    GQCP::AP1roGVariables variables (x, 2, 5);
    BOOST_CHECK(variables.asMatrix().isApprox(X));
}


BOOST_AUTO_TEST_CASE ( vectorIndex ) {

    size_t K = 11;
    size_t N_P = 2;
    GQCP::AP1roGVariables variables (N_P, K);

    BOOST_CHECK_EQUAL(variables.vectorIndex(0, 2), 0);
    BOOST_CHECK_EQUAL(variables.vectorIndex(0, 3), 1);
    BOOST_CHECK_EQUAL(variables.vectorIndex(1, 2), 9);
    BOOST_CHECK_EQUAL(variables.vectorIndex(1, 3), 10);

    // Require a throw if i > N_P
    BOOST_REQUIRE_THROW(variables.vectorIndex(3, 3), std::invalid_argument);

    // Require a throw if a < N_P
    BOOST_REQUIRE_THROW(variables.vectorIndex(0, 1), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( matrixIndex ) {

    size_t K = 11;
    size_t N_P = 2;
    GQCP::AP1roGVariables variables (N_P, K);


    BOOST_CHECK_EQUAL(variables.matrixIndexMajor(0), 0);
    BOOST_CHECK_EQUAL(variables.matrixIndexMajor(1), 0);
    BOOST_CHECK_EQUAL(variables.matrixIndexMajor(9), 1);
    BOOST_CHECK_EQUAL(variables.matrixIndexMajor(10), 1);

    BOOST_CHECK_EQUAL(variables.matrixIndexMinor(0), 2);
    BOOST_CHECK_EQUAL(variables.matrixIndexMinor(1), 3);
    BOOST_CHECK_EQUAL(variables.matrixIndexMinor(4), 6);
    BOOST_CHECK_EQUAL(variables.matrixIndexMinor(5), 7);
}


BOOST_AUTO_TEST_CASE ( operator_call ) {
    // Make an example for geminal coefficients for N_P=2 and K=5
    //      . .  1 2 3
    //      . .  4 5 6
    Eigen::VectorXd g (6);
    g << 1, 2, 3, 4, 5, 6;
    GQCP::AP1roGVariables variables (g, 2, 5);

    BOOST_CHECK_EQUAL(variables(0), 1);
    BOOST_CHECK_EQUAL(variables(1), 2);
    BOOST_CHECK_EQUAL(variables(2), 3);
    BOOST_CHECK_EQUAL(variables(3), 4);
    BOOST_CHECK_EQUAL(variables(4), 5);
    BOOST_CHECK_EQUAL(variables(5), 6);


    BOOST_CHECK_EQUAL(variables(0, 2), 1);
    BOOST_CHECK_EQUAL(variables(0, 3), 2);
    BOOST_CHECK_EQUAL(variables(0, 4), 3);
    BOOST_CHECK_EQUAL(variables(1, 2), 4);
    BOOST_CHECK_EQUAL(variables(1, 3), 5);
    BOOST_CHECK_EQUAL(variables(1, 4), 6);
}
