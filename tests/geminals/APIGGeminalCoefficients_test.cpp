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
#define BOOST_TEST_MODULE "APIGGeminalCoefficients"


#include "geminals/APIGGeminalCoefficients.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( numberOfGeminalCoefficients ) {

    BOOST_CHECK_EQUAL(GQCP::APIGGeminalCoefficients::numberOfGeminalCoefficients(2, 5), 10);

    BOOST_CHECK_THROW(GQCP::APIGGeminalCoefficients::numberOfGeminalCoefficients(4, 4), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( constructor ) {

    // Check a correct constructor
    GQCP::APIGGeminalCoefficients g (4, 6);

    // We can't create 4 geminals in 4 orbitals
    BOOST_CHECK_THROW(GQCP::APIGGeminalCoefficients (4, 4), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( constructor_vector ) {

    // Check a correct constructor
    Eigen::VectorXd g = Eigen::VectorXd::Zero(10);
    BOOST_CHECK_NO_THROW(GQCP::APIGGeminalCoefficients (g, 2, 5));

    // Check wrong parameters N_P and K
    BOOST_CHECK_THROW(GQCP::APIGGeminalCoefficients (g, 2, 4), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::APIGGeminalCoefficients (g, 1, 5), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( asMatrix ) {

    // For N_P=2 and K=5, we have an APIG geminal coefficient matrix that looks like the following matrix:
    Eigen::MatrixXd G (2, 5);
    G << 1, 2, 3, 4, 5,
         6, 7, 8, 9, 10;

    // The geminal coefficients, arranged in a vector, are then represented by the following vector:
    Eigen::VectorXd g (10);
    g << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;


    GQCP::APIGGeminalCoefficients gem_coeff (g, 2, 5);
    BOOST_CHECK(gem_coeff.asMatrix().isApprox(G));
}


BOOST_AUTO_TEST_CASE ( vectorIndex ) {

    size_t K = 5;
    size_t N_P = 2;
    GQCP::APIGGeminalCoefficients gem_coeff (N_P, K);


    BOOST_CHECK_EQUAL(gem_coeff.vectorIndex(0, 2), 2);
    BOOST_CHECK_EQUAL(gem_coeff.vectorIndex(1, 1), 6);
    BOOST_CHECK_EQUAL(gem_coeff.vectorIndex(1, 2), 7);

    // Require a throw if i > N_P
    BOOST_REQUIRE_THROW(gem_coeff.vectorIndex(3, 3), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( matrixIndex ) {

    size_t K = 5;
    size_t N_P = 2;
    GQCP::APIGGeminalCoefficients gem_coeff (N_P, K);


    BOOST_CHECK_EQUAL(gem_coeff.matrixIndexMajor(2), 0);
    BOOST_CHECK_EQUAL(gem_coeff.matrixIndexMajor(6), 1);
    BOOST_CHECK_EQUAL(gem_coeff.matrixIndexMajor(7), 1);

    BOOST_CHECK_EQUAL(gem_coeff.matrixIndexMinor(2), 2);
    BOOST_CHECK_EQUAL(gem_coeff.matrixIndexMinor(6), 1);
    BOOST_CHECK_EQUAL(gem_coeff.matrixIndexMinor(7), 2);
}


BOOST_AUTO_TEST_CASE ( operator_call ) {
    // Make an example for geminal coefficients for N_P=2 and K=5
    //      1 2 3 4 5
    //      6 7 8 9 10
    Eigen::VectorXd g (10);
    g << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
    GQCP::APIGGeminalCoefficients gem_coeff (g, 2, 5);


    BOOST_CHECK_EQUAL(gem_coeff(2), 3);
    BOOST_CHECK_EQUAL(gem_coeff(6), 7);
    BOOST_CHECK_EQUAL(gem_coeff(7), 8);

    BOOST_CHECK_EQUAL(gem_coeff(0, 2), 3);
    BOOST_CHECK_EQUAL(gem_coeff(1, 1), 7);
    BOOST_CHECK_EQUAL(gem_coeff(1, 2), 8);
}


BOOST_AUTO_TEST_CASE ( toWaveFunction_example1 ) {

    size_t K = 3;
    size_t N_P = 1;

    Eigen::VectorXd g (3);
    g << 1, 2, 3;
    GQCP::APIGGeminalCoefficients gem_coeff (g, N_P, K);

    Eigen::VectorXd ref_coefficients (3);
    ref_coefficients << 1, 2, 3;

BOOST_CHECK(ref_coefficients.isApprox(gem_coeff.toWaveFunction().get_coefficients()));
}


BOOST_AUTO_TEST_CASE ( toWaveFunction_example2 ) {

    size_t K = 3;
    size_t N_P = 2;

    Eigen::VectorXd g (6);
    g << 1, 2, 3, 4, 5, 6;
    GQCP::APIGGeminalCoefficients gem_coeff (g, N_P, K);

    Eigen::VectorXd ref_coefficients (3);
    ref_coefficients << 13, 18, 27;

    BOOST_CHECK(ref_coefficients.isApprox(gem_coeff.toWaveFunction().get_coefficients()));
}


BOOST_AUTO_TEST_CASE ( toWaveFunction_example3 ) {

    size_t K = 5;
    size_t N_P = 2;

    Eigen::VectorXd g (10);
    g << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
    GQCP::APIGGeminalCoefficients gem_coeff (g, N_P, K);

    Eigen::VectorXd ref_coefficients (10);
    ref_coefficients << 19, 26, 37, 33, 46, 59, 40, 55, 70, 85;

    BOOST_CHECK(ref_coefficients.isApprox(gem_coeff.toWaveFunction().get_coefficients()));
}
