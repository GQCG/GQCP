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
#define BOOST_TEST_MODULE "FrozenProductFockSpace"

#include <boost/test/unit_test.hpp>

#include "FockSpace/FrozenProductFockSpace.hpp"
#include "FockSpace/SelectedFockSpace.hpp"


BOOST_AUTO_TEST_CASE ( FrozenProductFockSpace_constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 1));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 2));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 3));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 4));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 5));

    GQCP::ProductFockSpace product_fock_space (10, 5, 5);

    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 1));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 2));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 3));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 4));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 5));
}


BOOST_AUTO_TEST_CASE ( FrozenProductFockSpace_member_test ) {

    GQCP::FrozenProductFockSpace frozen_space (10, 5, 5, 2);

    const GQCP::FrozenFockSpace& alpha_member = frozen_space.get_frozen_fock_space_alpha();
    const GQCP::FrozenFockSpace& beta_member = frozen_space.get_frozen_fock_space_beta();

    BOOST_CHECK(alpha_member.get_N() == 5);
    BOOST_CHECK(beta_member.get_N() == 5);

    BOOST_CHECK(alpha_member.get_K() == 10);
    BOOST_CHECK(beta_member.get_K() == 10);
}


BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_diagonal_vs_no_diagonal) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    auto parameters = GQCP::SQHamiltonian<double>::Molecular(hchain, "STO-3G");
    parameters.LowdinOrthonormalize();

    GQCP::FrozenProductFockSpace fock_space (6, 4, 4, 2);

    GQCP::SquareMatrix<double> hamiltonian = fock_space.evaluateOperatorDense(parameters, true);
    GQCP::SquareMatrix<double> hamiltonian_no_diagonal = fock_space.evaluateOperatorDense(parameters, false);
    GQCP::VectorX<double> hamiltonian_diagonal = fock_space.evaluateOperatorDiagonal(parameters);

    // Test if non-diagonal evaluation and diagonal evaluations are correct
    BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));
}


BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_Dense_diagonal_true ) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    auto parameters = GQCP::SQHamiltonian<double>::Molecular(hchain, "STO-3G");
    parameters.LowdinOrthonormalize();

    GQCP::FrozenProductFockSpace product_fock_space(6, 4, 4, 2);
    GQCP::SelectedFockSpace selected_fock_space(product_fock_space);

    auto& h = parameters.core();
    auto& g = parameters.twoElectron();

    // Test the evaluation of the operators with selected Fock space (the reference) versus the that of the product Fock space 
    auto one_electron_evaluation1 = product_fock_space.evaluateOperatorDense(h, true);
    auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(h, true);

    auto two_electron_evaluation1 = product_fock_space.evaluateOperatorDense(g, true);
    auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(g, true);

    auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDense(parameters, true);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(parameters, true);

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_Dense_diagonal_false ) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    auto parameters = GQCP::SQHamiltonian<double>::Molecular(hchain, "STO-3G");
    parameters.LowdinOrthonormalize();

    GQCP::FrozenProductFockSpace product_fock_space(6, 4, 4, 2);
    GQCP::SelectedFockSpace selected_fock_space(product_fock_space);

    auto& h = parameters.core();
    auto& g = parameters.twoElectron();

    auto one_electron_evaluation1 = product_fock_space.evaluateOperatorDense(h, false);
    auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(h, false);

    auto two_electron_evaluation1 = product_fock_space.evaluateOperatorDense(g, false);
    auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(g, false);

    auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDense(parameters, false);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(parameters, false);

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_diagonal ) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    auto parameters = GQCP::SQHamiltonian<double>::Molecular(hchain, "STO-3G");
    parameters.LowdinOrthonormalize();

    GQCP::FrozenProductFockSpace product_fock_space(6, 4, 4, 2);
    GQCP::SelectedFockSpace selected_fock_space(product_fock_space);

    auto& h = parameters.core();
    auto& g = parameters.twoElectron();

    auto one_electron_evaluation1 = product_fock_space.evaluateOperatorDiagonal(h);
    auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(h);

    auto two_electron_evaluation1 = product_fock_space.evaluateOperatorDiagonal(g);
    auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(g);

    auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDiagonal(parameters);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(parameters);

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}

