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

#define BOOST_TEST_MODULE "USQHamiltonian"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/USQHamiltonian.hpp"
#include "Utilities/linalg.hpp"
#include "Utilities/miscellaneous.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check the construction of USQHamiltonian, with faulty and correct inputs
 */
BOOST_AUTO_TEST_CASE(USQHamiltonian_constructor) {

    // Create single-particle basis
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::USpinorBasis<double, GQCP::GTOShell> spinor_basis {water, "STO-3G"};

    // Create One- and SQTwoElectronOperators (and a transformation matrix) with compatible dimensions
    const size_t K = spinor_basis.numberOfSpinors(GQCP::Spin::alpha);
    const GQCP::SquareMatrix<double> H_core = GQCP::SquareMatrix<double>::Random(K);
    GQCP::QCRankFourTensor<double> g {K};
    g.setRandom();

    // Create SQ operators with greater dimensions
    const GQCP::SquareMatrix<double> H_core_faulty = GQCP::SquareMatrix<double>::Random(K + 1);
    GQCP::QCRankFourTensor<double> g_faulty {K + 1};
    g_faulty.setRandom();

    // Create SQHamilonians with different dimensions
    const GQCP::SQHamiltonian<double> sq_hamiltonian_a {GQCP::ScalarSQOneElectronOperator<double>(H_core), GQCP::ScalarSQTwoElectronOperator<double>(g)};
    const GQCP::SQHamiltonian<double> sq_hamiltonian_b {GQCP::ScalarSQOneElectronOperator<double>(H_core), GQCP::ScalarSQTwoElectronOperator<double>(g)};
    const GQCP::SQHamiltonian<double> sq_hamiltonian_b_faulty {GQCP::ScalarSQOneElectronOperator<double>(H_core_faulty), GQCP::ScalarSQTwoElectronOperator<double>(g_faulty)};

    // Check if a correct constructor works with compatible elements
    BOOST_CHECK_NO_THROW(GQCP::USQHamiltonian<double> usq_hamiltonian(sq_hamiltonian_a, sq_hamiltonian_b, GQCP::ScalarSQTwoElectronOperator<double>(g)));
    // Check if a constructor throws an error with incompatible elements
    BOOST_CHECK_THROW(GQCP::USQHamiltonian<double> usq_hamiltonian(sq_hamiltonian_a, sq_hamiltonian_b_faulty, GQCP::ScalarSQTwoElectronOperator<double>(g)), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::USQHamiltonian<double> usq_hamiltonian(sq_hamiltonian_a, sq_hamiltonian_b, GQCP::ScalarSQTwoElectronOperator<double>(g_faulty)), std::invalid_argument);
}

/**
 *  Check if a total transformation or two individual transformations for the individual components of the USQHamiltonian amount to the same result
 */
BOOST_AUTO_TEST_CASE(USQHamiltonian_transform) {
    // Create single-particle basis for alpha and beta
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::USpinorBasis<double, GQCP::GTOShell> spinor_basis {water, "STO-3G"};

    const size_t K = spinor_basis.numberOfSpinors(GQCP::Spin::alpha);

    // Create two identical usq Hamiltonians
    auto usq_hamiltonian1 = GQCP::USQHamiltonian<double>::Molecular(spinor_basis, water);
    auto usq_hamiltonian2 = GQCP::USQHamiltonian<double>::Molecular(spinor_basis, water);

    const GQCP::SquareMatrix<double> U = GQCP::SquareMatrix<double>::RandomUnitary(K);

    // Perform a total transform and individual component transfromations
    usq_hamiltonian1.transform(U);
    usq_hamiltonian2.transform(U, GQCP::Spin::alpha);
    usq_hamiltonian2.transform(U, GQCP::Spin::beta);

    // Test if the transformation results in identical Hamtilonians
    BOOST_CHECK(usq_hamiltonian1.twoElectronMixed().parameters().isApprox(usq_hamiltonian2.twoElectronMixed().parameters()));
    BOOST_CHECK(usq_hamiltonian1.spinHamiltonian(GQCP::Spin::alpha).core().parameters().isApprox(usq_hamiltonian2.spinHamiltonian(GQCP::Spin::alpha).core().parameters()));
    BOOST_CHECK(usq_hamiltonian1.spinHamiltonian(GQCP::Spin::beta).core().parameters().isApprox(usq_hamiltonian2.spinHamiltonian(GQCP::Spin::beta).core().parameters()));
}
