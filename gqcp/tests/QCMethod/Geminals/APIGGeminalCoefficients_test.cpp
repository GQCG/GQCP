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
#define BOOST_TEST_MODULE "APIGGeminalCoefficients"

#include <boost/test/unit_test.hpp>

#include "QCMethod/Geminals/APIGGeminalCoefficients.hpp"


/**
 *  Check if the number of geminal coefficients is correctly implemented.
 */
BOOST_AUTO_TEST_CASE ( numberOfGeminalCoefficients ) {

    BOOST_CHECK_EQUAL(GQCP::APIGGeminalCoefficients::numberOfGeminalCoefficients(2, 5), 10);
    BOOST_CHECK_THROW(GQCP::APIGGeminalCoefficients::numberOfGeminalCoefficients(4, 4), std::invalid_argument);
}


/**
 *  Check if the construction of APIG geminal coefficients from a row-major vector represention is correct.
 */
BOOST_AUTO_TEST_CASE ( FromRowMajor ) {

    // For N_P=2 and K=5, we have an APIG geminal coefficient matrix that looks like the following matrix:
    GQCP::MatrixX<double> G (2, 5);
    G << 1, 2, 3, 4, 5,
         6, 7, 8, 9, 10;


    // The geminal coefficients, arranged in a vector, are then represented by the following vector:
    GQCP::VectorX<double> g (10);
    g << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;


    const auto gem_coeff = GQCP::APIGGeminalCoefficients::FromRowMajor(g, 2, 5);
    BOOST_CHECK(gem_coeff.asMatrix().isApprox(G));
}


/**
 *  Check if the construction of AP1roG geminal coefficients from a column-major vector represention is correct.
 */
BOOST_AUTO_TEST_CASE ( FromColumnMajor ) {

    // For N_P=2 and K=5, we have an AP1roG geminal coefficient matrix that looks like the following matrix:
    GQCP::MatrixX<double> G (2, 5);
    G << 1, 2, 3, 4, 5,
         6, 7, 8, 9, 10;


    // Test that we get the previous representation if we use the following vector that uses column-major indexing.
    GQCP::VectorX<double> g (10);
    g << 1, 6, 2, 7, 3, 8, 4, 9, 5, 10;


    const auto gem_coeff = GQCP::APIGGeminalCoefficients::FromColumnMajor(g, 2, 5);
    BOOST_CHECK(gem_coeff.asMatrix().isApprox(G));
}


/**
 *  Test if the conversion from APIG geminal coefficients to a wave function is correct (example 1).
 */
BOOST_AUTO_TEST_CASE ( toLinearExpansion_example1 ) {

    // Set up the normalized reference coefficients
    GQCP::VectorX<double> ref_coefficients (3);
    ref_coefficients << 1, 2, 3;
    ref_coefficients.normalize();


    // Construct the toy geminal coefficients
    const size_t K = 3;
    const size_t N_P = 1;

    GQCP::VectorX<double> g (3);
    g << 1, 2, 3;
    const auto gem_coeff = GQCP::APIGGeminalCoefficients::FromRowMajor(g, N_P, K);


    // Calculate the conversion from geminal coefficients to a wave function and check the result
    GQCP::SeniorityZeroONVBasis onv_basis (K, N_P);
    BOOST_CHECK(ref_coefficients.isApprox(gem_coeff.toLinearExpansion(onv_basis).coefficients()));
}


/**
 *  Test if the conversion from APIG geminal coefficients to a wave function is correct (example 2)
 */
BOOST_AUTO_TEST_CASE ( toLinearExpansion_example2 ) {

    // Set up the normalized reference coefficients
    GQCP::VectorX<double> ref_coefficients (3);
    ref_coefficients << 13, 18, 27;
    ref_coefficients.normalize();


    // Construct the toy geminal coefficients
    const size_t K = 3;
    const size_t N_P = 2;

    GQCP::VectorX<double> g (6);
    g << 1, 2, 3, 4, 5, 6;
    const auto gem_coeff = GQCP::APIGGeminalCoefficients::FromRowMajor(g, N_P, K);


    // Calculate the conversion from geminal coefficients to a wave function and check the result.
    GQCP::SeniorityZeroONVBasis onv_basis (K, N_P);
    BOOST_CHECK(ref_coefficients.isApprox(gem_coeff.toLinearExpansion(onv_basis).coefficients()));
}


/**
 *  Test if the conversion from APIG geminal coefficients to a wave function is correct (example 2).
 */
BOOST_AUTO_TEST_CASE ( toLinearExpansion_example3 ) {

    // Set up the normalized reference coefficients
    GQCP::VectorX<double> ref_coefficients (10);
    ref_coefficients << 19, 26, 37, 33, 46, 59, 40, 55, 70, 85;
    ref_coefficients.normalize();


    // Construct the toy geminal coefficients.
    const size_t K = 5;
    const size_t N_P = 2;

    GQCP::VectorX<double> g (10);
    g << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
    const auto gem_coeff = GQCP::APIGGeminalCoefficients::FromRowMajor(g, N_P, K);


    // Calculate the conversion from geminal coefficients to a wave function and check the result.
    GQCP::SeniorityZeroONVBasis onv_basis (K, N_P);
    BOOST_CHECK(ref_coefficients.isApprox(gem_coeff.toLinearExpansion(onv_basis).coefficients()));
}
