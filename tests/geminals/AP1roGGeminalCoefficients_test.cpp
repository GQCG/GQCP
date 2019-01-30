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
#define BOOST_TEST_MODULE "AP1roGGeminalCoefficients"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "geminals/AP1roGGeminalCoefficients.hpp"



BOOST_AUTO_TEST_CASE ( numberOfGeminalCoefficients ) {

    BOOST_CHECK_EQUAL(GQCP::AP1roGGeminalCoefficients::numberOfGeminalCoefficients(2, 5), 6);

    BOOST_CHECK_THROW(GQCP::AP1roGGeminalCoefficients::numberOfGeminalCoefficients(4, 4), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( asMatrix ) {

    // For N_P=2 and K=5, we have an AP1roG geminal coefficient matrix that looks like the following matrix:
    Eigen::MatrixXd G (2, 5);
    G << 1, 0,  1, 2, 3,
         0, 1,  4, 5, 6;

    // The geminal coefficients, arranged in a vector are then represented by the following vector:
    Eigen::VectorXd g (6);
    g << 1, 2, 3, 4, 5, 6;


    GQCP::AP1roGGeminalCoefficients gem_coeff (g, 2, 5);
    BOOST_CHECK(gem_coeff.asMatrix().isApprox(G));
}


BOOST_AUTO_TEST_CASE ( toWaveFunction_example1 ) {

    size_t K = 3;
    size_t N_P = 1;

    Eigen::VectorXd g (2);
    g << 2, 3;
    GQCP::AP1roGGeminalCoefficients gem_coeff (g, N_P, K);

    Eigen::VectorXd ref_coefficients (3);
    ref_coefficients << 1, 2, 3;

    GQCP::FockSpace fock_space (K, N_P);
    BOOST_CHECK(ref_coefficients.isApprox(gem_coeff.toWaveFunction(fock_space).get_coefficients()));
}


BOOST_AUTO_TEST_CASE ( toWaveFunction_example2 ) {

    size_t K = 3;
    size_t N_P = 2;

    Eigen::VectorXd g (2);
    g << 2, 3;
    GQCP::AP1roGGeminalCoefficients gem_coeff (g, N_P, K);

    Eigen::VectorXd ref_coefficients (3);
    ref_coefficients << 1, 3, 2;

    GQCP::FockSpace fock_space (K, N_P);
    BOOST_CHECK(ref_coefficients.isApprox(gem_coeff.toWaveFunction(fock_space).get_coefficients()));
}


BOOST_AUTO_TEST_CASE ( toWaveFunction_example3 ) {

    size_t K = 5;
    size_t N_P = 2;

    Eigen::VectorXd g (6);
    g << 2, 3, 4, 5, 6, 7;
    GQCP::AP1roGGeminalCoefficients gem_coeff (g, N_P, K);

    Eigen::VectorXd ref_coefficients (10);
    ref_coefficients << 1, 5, 2, 6, 3, 27, 7, 4, 34, 45;

    GQCP::FockSpace fock_space (K, N_P);
    BOOST_CHECK(ref_coefficients.isApprox(gem_coeff.toWaveFunction(fock_space).get_coefficients()));
}
