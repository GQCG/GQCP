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
#define BOOST_TEST_MODULE "GSpinorBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/GSpinorBasis.hpp"

#include "Molecule/Molecule.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"


/**
 *  Test if the most elementary constructor throws errors when expected
 */
BOOST_AUTO_TEST_CASE ( constructor_throws ) {

    // Initialize two scalar bases, one for the alpha component and one for the beta component
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");

    const GQCP::ScalarBasis<GQCP::GTOShell> alpha_scalar_basis (h2, "STO-3G");
    const auto K_alpha = alpha_scalar_basis.numberOfBasisFunctions();  // 2

    const GQCP::ScalarBasis<GQCP::GTOShell> beta_scalar_basis (h2, "6-31G");
    const auto K_beta = beta_scalar_basis.numberOfBasisFunctions();  // 4


    // Initialize two coefficient matrices, one with compatible and one with incompatible dimensions
    const GQCP::TransformationMatrix<double> T_compatible (K_alpha + K_beta);
    const GQCP::TransformationMatrix<double> T_incompatible (K_alpha + K_beta - 1);


    // Check if the constructor throws upon receiving incompatible arguments
    using SpinorBasisType = GQCP::GSpinorBasis<double, GQCP::GTOShell>;  // needed to resolve compilation errors with boost
    BOOST_CHECK_THROW(SpinorBasisType spinor_basis (alpha_scalar_basis, beta_scalar_basis, T_incompatible), std::invalid_argument);
    BOOST_CHECK_NO_THROW(SpinorBasisType spinor_basis (alpha_scalar_basis, beta_scalar_basis, T_compatible));
}


/**
 *  Test some basic functionality: number of spinors, number of alpha coefficients, number of beta coefficients
 */
BOOST_AUTO_TEST_CASE ( basic_functionality ) {

    // Initialize a spinor basis
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "STO-3G");


    // Check some basic functionality
    BOOST_CHECK(spinor_basis.numberOfAlphaCoefficients() == 2);
    BOOST_CHECK(spinor_basis.numberOfBetaCoefficients() == 2);
    BOOST_CHECK(spinor_basis.numberOfSpinors() == 4);
}


/**
 *  Check if the alpha- and beta-coefficient matrices together form the total coefficient matrix
 */
BOOST_AUTO_TEST_CASE ( alpha_beta_coefficient_matrix ) {

    // Initialize a spinor basis with a different scalar basis for both the components
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "STO-3G", "6-31G");
    const auto K_alpha = spinor_basis.numberOfAlphaCoefficients();
    const auto K_beta = spinor_basis.numberOfBetaCoefficients();
    const auto M = K_alpha + K_beta;  // number of spinors


    // Initialize reference values for the alpha, beta and total coefficient matrix and check the results
    GQCP::MatrixX<double> C_alpha_ref = GQCP::MatrixX<double>::Zero(K_alpha, M);
    C_alpha_ref.topLeftCorner(K_alpha, K_alpha) = GQCP::MatrixX<double>::Identity(K_alpha, K_alpha);

    GQCP::MatrixX<double> C_beta_ref = GQCP::MatrixX<double>::Zero(K_beta, M);
    C_beta_ref.bottomRightCorner(K_beta, K_beta) = GQCP::MatrixX<double>::Identity(K_beta, K_beta);

    GQCP::SquareMatrix<double> C_ref = GQCP::SquareMatrix<double>::Identity(M, M);

    BOOST_CHECK(spinor_basis.alphaCoefficientMatrix().isApprox(C_alpha_ref, 1.0e-08));
    BOOST_CHECK(spinor_basis.betaCoefficientMatrix().isApprox(C_beta_ref, 1.0e-08));
    BOOST_CHECK(spinor_basis.coefficientMatrix().isApprox(C_ref, 1.0e-08));
}
