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

#define BOOST_TEST_MODULE "SpinResolvedONVBasis"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/transform.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"


/**
 *  Test the SpinResolvedONVBasis constructor
 */
BOOST_AUTO_TEST_CASE(ProductONVBasis_constructor) {

    BOOST_CHECK_NO_THROW(GQCP::SpinResolvedONVBasis(10, 5, 5));
}


/**
 *  Check if the static SpinResolvedONVBasis basis dimension calculation is correct and if it can throw errors
 */
BOOST_AUTO_TEST_CASE(ProductONVBasis_dimension) {

    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(10, 1, 1), 100);
    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(6, 2, 2), 225);
    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(8, 3, 3), 3136);

    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(10, 2, 0), 45);
    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(6, 3, 1), 120);
    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(8, 4, 2), 1960);

    BOOST_CHECK_THROW(GQCP::SpinResolvedONVBasis::calculateDimension(60, 25, 25), std::overflow_error);
}


/**
 *  Perform a dense evaluation of a one-, two-electron operator and the Hamiltonian in the SpinResolvedONVBasis basis (including the diagonal)
 *  and compare these to the selected CI solutions.
 */
// BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_Dense_diagonal_true) {

//     GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
//     GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
//     spinor_basis.lowdinOrthonormalize();
//     auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

//     GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 4};
//     GQCP::SpinResolvedSelectedONVBasis selected_fock_space {product_fock_space};

//     const auto& h = sq_hamiltonian.core();
//     const auto& g = sq_hamiltonian.twoElectron();

//     // Test the evaluation of the operators with selected ONV basis (the reference) versus that of thhe specialized ONV basis
//     auto one_electron_evaluation1 = product_fock_space.evaluateOperatorDense(h, true);
//     auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(h, true);

//     auto two_electron_evaluation1 = product_fock_space.evaluateOperatorDense(g, true);
//     auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(g, true);

//     auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDense(sq_hamiltonian, true);
//     auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(sq_hamiltonian, true);

//     BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
//     BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
//     BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
// }


/**
 *  Perform a dense evaluation of a one-, two-electron operator and the Hamiltonian in the ONV basis (excluding the diagonal)
 *  and compare these to the selected CI solutions.
 */
// BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_Dense_diagonal_false) {

//     GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
//     GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
//     spinor_basis.lowdinOrthonormalize();
//     auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

//     GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 4};
//     GQCP::SpinResolvedSelectedONVBasis selected_fock_space {product_fock_space};

//     const auto& h = sq_hamiltonian.core();
//     const auto& g = sq_hamiltonian.twoElectron();

//     // Test the evaluation of the operators with selected ONV basis (the reference) versus that of the product ONV basis
//     auto one_electron_evaluation1 = product_fock_space.evaluateOperatorDense(h, false);
//     auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(h, false);

//     auto two_electron_evaluation1 = product_fock_space.evaluateOperatorDense(g, false);
//     auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(g, false);

//     auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDense(sq_hamiltonian, false);
//     auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(sq_hamiltonian, false);

//     BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
//     BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
//     BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
// }


/**
 *  Evaluate the diagonal of a one-, two-electron operator and the Hamiltonian in the ONV basis 
 *  and compare these to the selected CI solutions.
 */
// BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_diagonal) {

//     GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
//     GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
//     spinor_basis.lowdinOrthonormalize();
//     auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

//     GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 4};
//     GQCP::SpinResolvedSelectedONVBasis selected_fock_space {product_fock_space};

//     const auto& h = sq_hamiltonian.core();
//     const auto& g = sq_hamiltonian.twoElectron();

//     // Test the evaluation of the operators with selected ONV basis (the reference) versus that of the product ONV basis
//     auto one_electron_evaluation1 = product_fock_space.evaluateOperatorDiagonal(h);
//     auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(h);

//     auto two_electron_evaluation1 = product_fock_space.evaluateOperatorDiagonal(g);
//     auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(g);

//     auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDiagonal(sq_hamiltonian);
//     auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

//     BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
//     BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
//     BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
// }


/**
 *  Check the Dense evaluations with diagonal to that of the Dense with the diagonal excluded + the diagonal individually for the Hamiltonian
 */
// BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_diagonal_vs_no_diagonal) {

//     GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
//     GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
//     spinor_basis.lowdinOrthonormalize();
//     auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

//     GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 4};

//     GQCP::SquareMatrix<double> hamiltonian = product_fock_space.evaluateOperatorDense(sq_hamiltonian, true);
//     GQCP::SquareMatrix<double> hamiltonian_no_diagonal = product_fock_space.evaluateOperatorDense(sq_hamiltonian, false);
//     GQCP::VectorX<double> hamiltonian_diagonal = product_fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

//     // Test if non-diagonal evaluation and diagonal evaluations are correct
//     BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));
// }


/**
 *  Perform a matrix vector product evaluation of a one-, two-electron operator and the Hamiltonian in the ONV basis
 *  and compare these to the matrix vector product of the actual dense evaluations.
//  */
// BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_MatrixVectorProduct) {

//     GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
//     GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis(hchain, "STO-3G");
//     spinor_basis.lowdinOrthonormalize();
//     auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

//     GQCP::SpinResolvedONVBasis fock_space {6, 4, 4};

//     const auto& h = sq_hamiltonian.core();
//     const auto& g = sq_hamiltonian.twoElectron();

//     // Generate diagonals for the matvec input
//     auto one_electron_diagonal = fock_space.evaluateOperatorDiagonal(h);
//     auto two_electron_diagonal = fock_space.evaluateOperatorDiagonal(g);
//     auto hamiltonian_diagonal = fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

//     // Test the evaluation of the operators with selected ONV basis (the reference) versus that of the product ONV basis
//     auto one_electron_evaluation1 = fock_space.evaluateOperatorMatrixVectorProduct(h, one_electron_diagonal, one_electron_diagonal);
//     GQCP::VectorX<double> one_electron_evaluation2 = fock_space.evaluateOperatorDense(h, true) * one_electron_diagonal;

//     auto two_electron_evaluation1 = fock_space.evaluateOperatorMatrixVectorProduct(g, two_electron_diagonal, two_electron_diagonal);
//     GQCP::VectorX<double> two_electron_evaluation2 = fock_space.evaluateOperatorDense(g, true) * two_electron_diagonal;

//     auto hamiltonian_evaluation1 = fock_space.evaluateOperatorMatrixVectorProduct(sq_hamiltonian, hamiltonian_diagonal, hamiltonian_diagonal);
//     GQCP::VectorX<double> hamiltonian_evaluation2 = fock_space.evaluateOperatorDense(sq_hamiltonian, true) * hamiltonian_diagonal;

//     BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
//     BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
//     BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
// }


/**
 *  This tests the results for diagonal and dense evaluations for the restricted framework to that of the unrestricted framework in a restricted basis (the alpha and beta coefficients and parameters are identical) 
 */
// BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_diagonal_unrestricted) {

//     GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);

//     GQCP::USpinorBasis<double, GQCP::GTOShell> uspinor_basis {hchain, "STO-3G"};
//     GQCP::RSpinorBasis<double, GQCP::GTOShell> rspinor_basis {hchain, "STO-3G"};
//     uspinor_basis.lowdinOrthonormalize();
//     rspinor_basis.lowdinOrthonormalize();

//     auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(uspinor_basis, hchain);  // unrestricted Hamiltonian in the Löwdin basis
//     auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(rspinor_basis, hchain);   // restricted Hamiltonian in the Löwdin basis

//     GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 4};

//     auto hamiltonian_diagonal_evaluation1 = product_fock_space.evaluateOperatorDiagonal(sq_hamiltonian);
//     auto hamiltonian_diagonal_evaluation2 = product_fock_space.evaluateOperatorDiagonal(usq_hamiltonian);

//     auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDense(sq_hamiltonian, false);
//     auto hamiltonian_evaluation2 = product_fock_space.evaluateOperatorDense(usq_hamiltonian, false);

//     BOOST_CHECK(hamiltonian_diagonal_evaluation1.isApprox(hamiltonian_diagonal_evaluation2));
//     BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
// }


/**
 *  Perform a dense and diagonal evaluation for the unrestricted Hamiltonian in the product ONV basis and compare these to the selected CI solutions
 */
// BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_diagonal_unrestricted_vs_selected) {

//     GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
//     GQCP::USpinorBasis<double, GQCP::GTOShell> uspinor_basis {hchain, "STO-3G"};
//     uspinor_basis.lowdinOrthonormalize();
//     auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(uspinor_basis, hchain);  // restricted Hamiltonian in the Löwdin basis

//     // Transform the beta component
//     // Create stable unitairy matrix
//     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes {usq_hamiltonian.spinHamiltonian(GQCP::Spin::alpha).core().parameters()};
//     GQCP::basisTransform(uspinor_basis, usq_hamiltonian, GQCP::TransformationMatrix<double>(saes.eigenvectors()), GQCP::Spin::beta);

//     GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 4};
//     GQCP::SpinResolvedSelectedONVBasis selected_fock_space {product_fock_space};

//     auto hamiltonian_diagonal_evaluation1 = product_fock_space.evaluateOperatorDiagonal(usq_hamiltonian);
//     auto hamiltonian_diagonal_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(usq_hamiltonian);

//     auto hamiltonian_evaluation1 = product_fock_space.evaluateOperatorDense(usq_hamiltonian, true);
//     auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(usq_hamiltonian, true);

//     BOOST_CHECK(hamiltonian_diagonal_evaluation1.isApprox(hamiltonian_diagonal_evaluation2));
//     BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
// }
