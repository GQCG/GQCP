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

#include "Basis/transform.hpp"
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


/**
 *  Perform a dense evaluation of a one-, two-electron operator and the Hamiltonian in the frozen product Fock space (excluding the diagonal)
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_diagonal_vs_no_diagonal) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::FrozenProductFockSpace fock_space (6, 4, 4, 2);

    GQCP::SquareMatrix<double> hamiltonian = fock_space.evaluateOperatorDense(sq_hamiltonian, true);
    GQCP::SquareMatrix<double> hamiltonian_no_diagonal = fock_space.evaluateOperatorDense(sq_hamiltonian, false);
    GQCP::VectorX<double> hamiltonian_diagonal = fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

    // Test if non-diagonal evaluation and diagonal evaluations are correct
    BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));
}


/**
 *  Perform a dense evaluation of a one-, two-electron operator and the Hamiltonian in the frozen product Fock space (including the diagonal)
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_Dense_diagonal_true ) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::FrozenProductFockSpace product_fock_space(6, 4, 4, 2);
    GQCP::SelectedFockSpace selected_fock_space(product_fock_space);

    auto& h = sq_hamiltonian.core();
    auto& g = sq_hamiltonian.twoElectron();

    // Test the evaluation of the operators with selected Fock space (the reference) versus the that of the product Fock space 
    auto one_electron_evaluation1 = product_fock_space.evaluateOperatorDense(h, true);
    auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(h, true);

    auto two_electron_evaluation1 = product_fock_space.evaluateOperatorDense(g, true);
    auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(g, true);

    auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDense(sq_hamiltonian, true);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(sq_hamiltonian, true);

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


/**
 *  Perform a dense evaluation of a one-, two-electron operator and the Hamiltonian in the frozen product Fock space (excluding the diagonal)
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_Dense_diagonal_false ) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::FrozenProductFockSpace product_fock_space(6, 4, 4, 2);
    GQCP::SelectedFockSpace selected_fock_space(product_fock_space);

    auto& h = sq_hamiltonian.core();
    auto& g = sq_hamiltonian.twoElectron();

    auto one_electron_evaluation1 = product_fock_space.evaluateOperatorDense(h, false);
    auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(h, false);

    auto two_electron_evaluation1 = product_fock_space.evaluateOperatorDense(g, false);
    auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(g, false);

    auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDense(sq_hamiltonian, false);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(sq_hamiltonian, false);

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


/**
 *  Perform a diagonal evaluation of a one-, two-electron operator and the Hamiltonian in the frozen product Fock space
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_diagonal ) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::FrozenProductFockSpace product_fock_space(6, 4, 4, 2);
    GQCP::SelectedFockSpace selected_fock_space(product_fock_space);

    auto& h = sq_hamiltonian.core();
    auto& g = sq_hamiltonian.twoElectron();

    auto one_electron_evaluation1 = product_fock_space.evaluateOperatorDiagonal(h);
    auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(h);

    auto two_electron_evaluation1 = product_fock_space.evaluateOperatorDiagonal(g);
    auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(g);

    auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDiagonal(sq_hamiltonian);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


/**
 *  Perform a dense and diagonal evaluation for the unrestricted Hamiltonian in the frozen product Fock space
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE ( FrozenProductFockSpace_evaluateOperator_diagonal_unrestricted_vs_selected ) {

    // This test the evaluations in an unrestricted basis for the FrozenProductFockSpace versus the evaluation of the selected module
    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);

    GQCP::USpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    spinor_basis.lowdinOrthonormalize();

    auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(spinor_basis, hchain);  // unrestricted Hamiltonian in the Löwdin basis

    // Transform the beta component
    // Create stable unitairy matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (usq_hamiltonian.spinHamiltonian(GQCP::SpinComponent::ALPHA).core().parameters());
    GQCP::basisTransform(spinor_basis, usq_hamiltonian, GQCP::TransformationMatrix<double>(saes.eigenvectors()), GQCP::SpinComponent::BETA);

    GQCP::FrozenProductFockSpace product_fock_space (6, 4, 4, 2);
    GQCP::SelectedFockSpace selected_fock_space (product_fock_space);

    auto hamiltonian_diagonal_evaluation1 = product_fock_space.evaluateOperatorDiagonal(usq_hamiltonian);
    auto hamiltonian_diagonal_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(usq_hamiltonian);

    auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDense(usq_hamiltonian, true);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(usq_hamiltonian, true);

    BOOST_CHECK(hamiltonian_diagonal_evaluation1.isApprox(hamiltonian_diagonal_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}