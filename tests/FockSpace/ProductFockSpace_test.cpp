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
#define BOOST_TEST_MODULE "ProductFockSpace"

#include <boost/test/unit_test.hpp>

#include "FockSpace/ProductFockSpace.hpp"
#include "FockSpace/SelectedFockSpace.hpp"


BOOST_AUTO_TEST_CASE ( ProductFockSpace_constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::ProductFockSpace (10, 5, 5));
}


BOOST_AUTO_TEST_CASE ( ProductFockSpace_dimension) {

    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(10, 1, 1), 100);
    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(6, 2, 2), 225);
    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(8, 3, 3), 3136);

    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(10, 2, 0), 45);
    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(6, 3, 1), 120);
    BOOST_CHECK_EQUAL(GQCP::ProductFockSpace::calculateDimension(8, 4, 2), 1960);

    BOOST_CHECK_THROW(GQCP::ProductFockSpace::calculateDimension(60, 25, 25), std::overflow_error);

}


BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_Dense_diagonal_true ) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (hchain, "STO-3G");
    sp_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, hchain);  // in the Löwdin basis

    GQCP::ProductFockSpace product_fock_space (6, 4, 4);
    GQCP::SelectedFockSpace selected_fock_space (product_fock_space);

    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    // Test the evaluation of the operators with selected Fock space (the reference) versus that of the product Fock space 
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


BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_Dense_diagonal_false ) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (hchain, "STO-3G");
    sp_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, hchain);  // in the Löwdin basis

    GQCP::ProductFockSpace product_fock_space (6, 4, 4);
    GQCP::SelectedFockSpace selected_fock_space (product_fock_space);

    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    // Test the evaluation of the operators with selected Fock space (the reference) versus that of the product Fock space 
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


BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_diagonal ) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (hchain, "STO-3G");
    sp_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, hchain);  // in the Löwdin basis

    GQCP::ProductFockSpace product_fock_space (6, 4, 4);
    GQCP::SelectedFockSpace selected_fock_space (product_fock_space);

    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    // Test the evaluation of the operators with selected Fock space (the reference) versus that of the product Fock space 
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


BOOST_AUTO_TEST_CASE ( FockSpace_EvaluateOperator_diagonal_vs_no_diagonal) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (hchain, "STO-3G");
    sp_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, hchain);  // in the Löwdin basis

    GQCP::ProductFockSpace product_fock_space (6, 4, 4);

    GQCP::SquareMatrix<double> hamiltonian = product_fock_space.evaluateOperatorDense(sq_hamiltonian, true);
    GQCP::SquareMatrix<double> hamiltonian_no_diagonal = product_fock_space.evaluateOperatorDense(sq_hamiltonian, false);
    GQCP::VectorX<double> hamiltonian_diagonal = product_fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

    // Test if non-diagonal evaluation and diagonal evaluations are correct
    BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));
}
