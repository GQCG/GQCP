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

#define BOOST_TEST_MODULE "USpinOrbitalBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Molecule/Molecule.hpp"


/**
 *  Test if the most elementary constructor throws errors when expected.
 */
BOOST_AUTO_TEST_CASE(constructor_throws) {

    // Initialize two scalar bases, one for the alpha component and one for the beta component.
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");

    const GQCP::ScalarBasis<GQCP::GTOShell> alpha_scalar_basis {h2, "STO-3G"};
    const auto K_alpha = alpha_scalar_basis.numberOfBasisFunctions();  // 2

    const GQCP::ScalarBasis<GQCP::GTOShell> beta_scalar_basis {h2, "6-31G"};
    const auto K_beta = beta_scalar_basis.numberOfBasisFunctions();  // 4


    // Initialize compatible and incompatible test expansion coefficients.
    const GQCP::UTransformation<double> T_compatible {GQCP::UTransformationComponent<double>::Zero(K_alpha), GQCP::UTransformationComponent<double>::Zero(K_beta)};

    const GQCP::UTransformation<double> T_incompatible_alpha {GQCP::UTransformationComponent<double>::Zero(K_alpha + 1), GQCP::UTransformationComponent<double>::Zero(K_beta)};
    const GQCP::UTransformation<double> T_incompatible_beta {GQCP::UTransformationComponent<double>::Zero(K_alpha), GQCP::UTransformationComponent<double>::Zero(K_beta + 1)};


    // Check if the constructor throws upon receiving incompatible arguments.
    using SpinorBasisType = GQCP::USpinOrbitalBasis<double, GQCP::GTOShell>;  // needed to resolve compilation errors with boost
    BOOST_CHECK_THROW(SpinorBasisType spinor_basis(alpha_scalar_basis, beta_scalar_basis, T_incompatible_alpha), std::invalid_argument);
    BOOST_CHECK_THROW(SpinorBasisType spinor_basis(alpha_scalar_basis, beta_scalar_basis, T_incompatible_beta), std::invalid_argument);

    BOOST_CHECK_NO_THROW(SpinorBasisType spinor_basis(alpha_scalar_basis, beta_scalar_basis, T_compatible));
}


/**
 *  Test some basic functionality: number of spinors/spin-orbitals.
 */
BOOST_AUTO_TEST_CASE(basic_functionality) {

    // Initialize an unrestricted spin-orbital basis.
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {h2, "STO-3G"};

    // Check if the number of spin-orbitals is correct.
    BOOST_CHECK(spinor_basis.numberOfSpinors() == 4);
}
