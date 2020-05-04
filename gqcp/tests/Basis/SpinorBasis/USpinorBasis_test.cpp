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

#define BOOST_TEST_MODULE "USpinorBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/USpinorBasis.hpp"
#include "Molecule/Molecule.hpp"


/**
 *  Test if the most elementary constructor throws errors when expected
 */
BOOST_AUTO_TEST_CASE(constructor_throws) {

    // Initialize two scalar bases, one for the alpha component and one for the beta component
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");

    const GQCP::ScalarBasis<GQCP::GTOShell> alpha_scalar_basis {h2, "STO-3G"};
    const auto K_alpha = alpha_scalar_basis.numberOfBasisFunctions();  // 2

    const GQCP::ScalarBasis<GQCP::GTOShell> beta_scalar_basis {h2, "6-31G"};
    const auto K_beta = beta_scalar_basis.numberOfBasisFunctions();  // 4


    // Initialize four coefficient matrices, two for each component, of which one with compatible and one with incompatible dimensions
    const GQCP::TransformationMatrix<double> T_compatible_alpha {K_alpha};
    const GQCP::TransformationMatrix<double> T_compatible_beta {K_beta};
    const GQCP::TransformationMatrix<double> T_incompatible_alpha {K_alpha - 1};
    const GQCP::TransformationMatrix<double> T_incompatible_beta {K_beta - 1};


    // Check if the constructor throws upon receiving incompatible arguments
    using SpinorBasisType = GQCP::USpinorBasis<double, GQCP::GTOShell>;  // needed to resolve compilation errors with boost
    BOOST_CHECK_THROW(SpinorBasisType spinor_basis(alpha_scalar_basis, beta_scalar_basis, T_compatible_alpha, T_incompatible_beta), std::invalid_argument);
    BOOST_CHECK_THROW(SpinorBasisType spinor_basis(alpha_scalar_basis, beta_scalar_basis, T_incompatible_alpha, T_compatible_beta), std::invalid_argument);
    BOOST_CHECK_THROW(SpinorBasisType spinor_basis(alpha_scalar_basis, beta_scalar_basis, T_incompatible_alpha, T_incompatible_beta), std::invalid_argument);
    BOOST_CHECK_NO_THROW(SpinorBasisType spinor_basis(alpha_scalar_basis, beta_scalar_basis, T_compatible_alpha, T_compatible_beta));
}


/**
 *  Test some basic functionality: number of spinors, number of alpha coefficients, number of beta coefficients
 */
BOOST_AUTO_TEST_CASE(basic_functionality) {

    // Initialize a spinor basis
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const GQCP::USpinorBasis<double, GQCP::GTOShell> spinor_basis {h2, "STO-3G"};


    // Check some basic functionality
    BOOST_CHECK(spinor_basis.numberOfCoefficients(GQCP::Spin::alpha) == 2);
    BOOST_CHECK(spinor_basis.numberOfCoefficients(GQCP::Spin::beta) == 2);
    BOOST_CHECK(spinor_basis.numberOfSpinors() == 4);
}
