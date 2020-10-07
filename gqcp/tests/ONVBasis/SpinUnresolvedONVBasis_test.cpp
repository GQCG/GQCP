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

#define BOOST_TEST_MODULE "SpinUnresolvedONVBasis"

#include <boost/test/unit_test.hpp>

#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


/**
 *  Check if the calculation of the dimension of a SpinUnresolvedONVBasis is correct and if it can throw errors.
 */
BOOST_AUTO_TEST_CASE(ONVBasis_dimension) {

    BOOST_CHECK_EQUAL(GQCP::SpinUnresolvedONVBasis::calculateDimension(10, 1), 10);
    BOOST_CHECK_EQUAL(GQCP::SpinUnresolvedONVBasis::calculateDimension(6, 2), 15);
    BOOST_CHECK_EQUAL(GQCP::SpinUnresolvedONVBasis::calculateDimension(8, 3), 56);

    BOOST_CHECK_THROW(GQCP::SpinUnresolvedONVBasis::calculateDimension(100, 50), std::overflow_error);
}


/**
 *  Test if the vertex weights of the SpinUnresolvedONV basis addressing scheme are correct for the Fock space F(5,3).
 * 
 *  For this (full) Fock space, the vertex weights are as follows:
 * 
 *      1   0   0   0
 *      1   1   0   0
 *      1   2   1   0
 *      0   3   3   1
 *      0   0   6   4
 *      0   0   0   10
 */
BOOST_AUTO_TEST_CASE(vertex_weights_M5_N3) {

    const GQCP::SpinUnresolvedONVBasis onv_basis {5, 3};

    BOOST_CHECK(onv_basis.vertexWeight(0, 0) == 1);
    BOOST_CHECK(onv_basis.vertexWeight(1, 2) == 0);
    BOOST_CHECK(onv_basis.vertexWeight(3, 2) == 3);
    BOOST_CHECK(onv_basis.vertexWeight(4, 2) == 6);
    BOOST_CHECK(onv_basis.vertexWeight(5, 3) == 10);
}


/**
 *  Test if the arc weights of the SpinUnresolvedONV basis addressing scheme are correct for the Fock space F(5,3).
 * 
 *  For this (full) Fock space, the vertex weights are as follows:
 * 
 *      1   0   0   0
 *      1   1   0   0
 *      1   2   1   0
 *      0   3   3   1
 *      0   0   6   4
 *      0   0   0   10
 */
BOOST_AUTO_TEST_CASE(arc_weights_M5_N3) {

    const GQCP::SpinUnresolvedONVBasis onv_basis {5, 3};

    BOOST_CHECK(onv_basis.arcWeight(3, 1) == 3);
    BOOST_CHECK(onv_basis.arcWeight(4, 2) == 4);
    BOOST_CHECK(onv_basis.arcWeight(2, 2) == 0);
}


/**
 *  Test if the shift in address and orbital indices correspond to the correct solutions
 *  
 *  This method does not alter the ONVs in any way.
 */
BOOST_AUTO_TEST_CASE(iterateToNextUnoccupiedOrbital) {

    GQCP::SpinUnresolvedONVBasis fock_space {5, 3};
    GQCP::SpinUnresolvedONV onv = fock_space.constructONVFromAddress(3);  // 01110

    size_t address_shift = 0;
    // test shift if we annihilate one electron and start from orbital index 2
    size_t e = 1;  // count starts at 1 (translates to orbital index 2)
    size_t q = 2;  // index starts at orbital index 2

    //  In this instance electron weights at index 2 and 3 should be shifted.
    //  Initial weight contributions were 1 and 1 respectively,
    //  these should be shifted to 2 and 3 respectively, the difference is 1 and 2 respectively.
    //  The total shift is thus 3
    fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address_shift, q, e);

    BOOST_CHECK(address_shift == 3);
    BOOST_CHECK(e == 3);
    BOOST_CHECK(q == 4);

    // test shift if we annihilate two electrons and start from orbital index 3
    e = 2;  // count starts at 2 (translates to orbital index 3)
    q = 3;  // index starts at orbital index 3

    //  In this instance electron weights at index 3 should be shifted.
    //  The initial weight contribution was 1,
    //  this should be shifted to 3, the difference is 2
    //  The total shift is thus 2
    address_shift = 0;
    fock_space.shiftUntilNextUnoccupiedOrbital<2>(onv, address_shift, q, e);

    BOOST_CHECK(address_shift == 2);
    BOOST_CHECK(e == 3);
    BOOST_CHECK(q == 4);
}


/**
 *  Test if the shift in address, orbital indices and sign change correspond to the correct solutions
 *  
 *  This method does not alter the ONVs in any way.
 */
BOOST_AUTO_TEST_CASE(iterateToNextUnoccupiedOrbital_signed) {

    GQCP::SpinUnresolvedONVBasis fock_space {5, 3};
    GQCP::SpinUnresolvedONV onv = fock_space.constructONVFromAddress(3);  // 01110

    size_t address_shift = 0;
    int sign = 1;
    // test shift if we annihilate one electron and start from orbital index 2
    size_t e = 1;  // count starts at 1 (translates to orbital index 2)
    size_t q = 2;  // index starts at orbital index 2

    //  In this instance electron weights at index 2 and 3 should be shifted.
    //  Initial weight contributions were 1 and 1 respectively,
    //  these should be shifted to 2 and 3 respectively, the difference is 1 and 2 respectively.
    //  The total shift is thus 3, and the sign should remain the same (flips twice)
    fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == 3);
    BOOST_CHECK(e == 3);
    BOOST_CHECK(q == 4);
    BOOST_CHECK(sign == 1);

    // test shift if we annihilate two electrons and start from orbital index 3
    e = 2;  // count starts at 2 (translates to orbital index 3)
    q = 3;  // index starts at orbital index 3

    //  In this instance electron weights at index 3 should be shifted.
    //  The initial weight contribution was 1,
    //  this should be shifted to 3, the difference is 2
    //  The total shift is thus 2, sign should flip once
    address_shift = 0;
    fock_space.shiftUntilNextUnoccupiedOrbital<2>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == 2);
    BOOST_CHECK(e == 3);
    BOOST_CHECK(q == 4);
    BOOST_CHECK(sign == -1);
}


/**
 *  Test if the shift in address, orbital indices and sign change correspond to the correct solutions
 *  
 *  This method does not alter the ONVs in any way.
 */
BOOST_AUTO_TEST_CASE(shiftToPreviousOrbital_signed) {

    GQCP::SpinUnresolvedONVBasis fock_space {5, 3};
    GQCP::SpinUnresolvedONV onv = fock_space.constructONVFromAddress(6);  // 10110

    size_t address_shift = 0;
    int sign = 1;

    // test shift if we plan on creating one electron and start from orbital index 2
    size_t e = 0;  // count starts at 0 (translates to orbital index 1)
    size_t q = 1;  // index starts at orbital index 1

    //  Index 1 is occupied an thus its weight shall shift if we create before an electron on a smaller index.
    //  In this instance the electron weight at index 1 should be shifted.
    //  Initial weight contribution is 1,
    //  This should be shifted to 0, the difference is 1.
    //  The sign changes (flips once)
    fock_space.shiftUntilPreviousUnoccupiedOrbital<1>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == -1);
    BOOST_CHECK(e == -1);
    BOOST_CHECK(q == 0);
    BOOST_CHECK(sign == -1);

    sign = 1;
    address_shift = 0;
    onv = fock_space.constructONVFromAddress(9);  // 11100
    // test shift if we plan on creating two electrons and start from orbital index 2
    e = 0;  // count starts at 1 (translates to orbital index 2)
    q = 2;  // index starts at orbital index 2

    //  Index 2 is occupied an thus its weight shall shift if we create before an electron on a smaller index.
    //  In this instance the electron weight at index 2 should be shifted.
    //  Initial weight contribution is 2,
    //  This should be shifted to 0, the difference is 2.
    //  The sign changes (flips once)
    fock_space.shiftUntilPreviousUnoccupiedOrbital<2>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == -2);
    BOOST_CHECK(e == -1);
    BOOST_CHECK(q == 1);
    BOOST_CHECK(sign == -1);
}


/**
 *  Test if the SpinUnresolvedONV basis correctly counts the number of coupling ONVs with larger address for a given SpinUnresolvedONV
 */
BOOST_AUTO_TEST_CASE(coupling_count) {

    GQCP::SpinUnresolvedONVBasis fock_space {5, 3};
    GQCP::SpinUnresolvedONV onv = fock_space.constructONVFromAddress(3);  // 01110

    // We only count couplings with larger addresses

    BOOST_CHECK(fock_space.countOneElectronCouplings(onv) == 3);      // 11100, 11010, 10110
    BOOST_CHECK(fock_space.countTwoElectronCouplings(onv) == 3 + 3);  // 11100, 11010, 10110, 11001, 10101, 10011

    onv = fock_space.constructONVFromAddress(0);  // 00111

    BOOST_CHECK(fock_space.countOneElectronCouplings(onv) == 6);
    BOOST_CHECK(fock_space.countTwoElectronCouplings(onv) == 6 + 3);  // all of them


    // test whether the total count matches that of individual counts of all ONVs in the SpinUnresolvedONV basis.
    GQCP::SpinUnresolvedONVBasis fock_space2 {16, 8};

    size_t coupling_count1 = 0;
    size_t coupling_count2 = 0;
    onv = fock_space2.constructONVFromAddress(0);           // spin string with address 0
    for (size_t I = 0; I < fock_space2.dimension(); I++) {  // I_alpha loops over all addresses of alpha spin strings
        if (I > 0) {
            fock_space2.transformONVToNextPermutation(onv);
        }
        coupling_count1 += fock_space2.countOneElectronCouplings(onv);
        coupling_count2 += fock_space2.countTwoElectronCouplings(onv);
    }

    BOOST_CHECK(2 * coupling_count1 == fock_space2.countTotalOneElectronCouplings());
    BOOST_CHECK(2 * coupling_count2 == fock_space2.countTotalTwoElectronCouplings());
}


/**
 *  In this test we iterate over the entire SpinUnresolvedONV basis using the SpinUnresolvedONVBasis::setNextONV(SpinUnresolvedONV&) 
 *  and test wether the address is correct using SpinUnresolvedONVBasis::address(const SpinUnresolvedONV&) 
 */
BOOST_AUTO_TEST_CASE(ONV_address_setNext_fullspace) {

    // Here we will test a full permutation through a SpinUnresolvedONV basis of K = 15, N = 5
    GQCP::SpinUnresolvedONVBasis fock_space {15, 5};

    // Retrieve the first SpinUnresolvedONV of the SpinUnresolvedONV basis
    GQCP::SpinUnresolvedONV onv_test = fock_space.constructONVFromAddress(0);

    const size_t dimension_fock_space = 3003;
    bool is_correct = true;  // variable that is updated to false if an unexpected result occurs

    // Iterate through the SpinUnresolvedONV basis in reverse lexicographical order and test whether address matches
    for (size_t i = 0; i < dimension_fock_space; i++) {

        // Tests address
        if (i != fock_space.addressOf(onv_test)) {
            is_correct = false;
        }

        // transforms the given SpinUnresolvedONV to the next SpinUnresolvedONV in the SpinUnresolvedONV basis
        if (i < dimension_fock_space - 1) {
            fock_space.transformONVToNextPermutation(onv_test);
        }
    }

    // Checks if no unexpected results occured in a full iteration
    BOOST_CHECK(is_correct);
}


/**
 *  Test wether the SpinUnresolvedONV basis attributes the correct address to a given SpinUnresolvedONV.
 */
BOOST_AUTO_TEST_CASE(ONVBasis_getAddress) {

    GQCP::SpinUnresolvedONVBasis fock_space {6, 3};

    // The address of the string "010011" (19) should be 4
    GQCP::SpinUnresolvedONV onv {6, 3, 19};

    BOOST_CHECK_EQUAL(fock_space.addressOf(onv), 4);
}


/**
 *  Test setNext for manually chosen ONVs
 */
BOOST_AUTO_TEST_CASE(ONVBasis_setNext) {

    GQCP::SpinUnresolvedONVBasis fock_space {5, 3};
    // K = 5, N = 3 <-> "00111"
    GQCP::SpinUnresolvedONV onv = fock_space.constructONVFromAddress(0);
    // The lexical permutations are: "00111" (7), "01011" (11), "01101" (13), "01110" (14), etc.

    // Check permutations one after the other

    fock_space.transformONVToNextPermutation(onv);  // "01011" (11)
    BOOST_CHECK_EQUAL(onv.unsignedRepresentation(), 11);
    std::vector<size_t> ref_indices1 {0, 1, 3};
    BOOST_CHECK(ref_indices1 == onv.occupiedIndices());

    fock_space.transformONVToNextPermutation(onv);  // "01101" (13)
    BOOST_CHECK_EQUAL(onv.unsignedRepresentation(), 13);
    std::vector<size_t> ref_indices2 {0, 2, 3};
    BOOST_CHECK(ref_indices2 == onv.occupiedIndices());

    fock_space.transformONVToNextPermutation(onv);  // "01110" (14)
    BOOST_CHECK_EQUAL(onv.unsignedRepresentation(), 14);
    std::vector<size_t> ref_indices3 {1, 2, 3};
    BOOST_CHECK(ref_indices3 == onv.occupiedIndices());
}


/**
 *  Check if the evaluation of a one-electron operator in a spin-unresolved ONV basis matches the evaluation in a spin-resolved selected ONV basis with only alpha electrons.
 * 
 *  The test system is a H6(2+)-chain with internuclear separation of 0.742 (a.u.) in an STO-3G basis.
 */
BOOST_AUTO_TEST_CASE(evaluate_one_electron_operator) {

    // Set up an example molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HChain(6, 0.742, +2);
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto M = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in the orthonormal Löwdin basis

    // Set up the two equivalent ONV bases.
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};
    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {GQCP::SpinResolvedONVBasis(M, N, 0)};  // no beta electrons


    // Check the evaluation of the core Hamiltonian.
    const auto& h_op = sq_hamiltonian.core();

    // Check the dense evaluation.
    const auto h_dense = onv_basis.evaluate<GQCP::SquareMatrix<double>>(h_op, true);     // true: calculate diagonal values
    const auto h_dense_selected = selected_onv_basis.evaluateOperatorDense(h_op, true);  // true: calculate diagonal values

    BOOST_CHECK(h_dense.isApprox(h_dense_selected));


    // Check the sparse evaluation.


    //
}


/**
 *  Perform a dense evaluation of a one-, two-electron operator and the Hamiltonian in the SpinUnresolvedONV basis (including the diagonal)
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_Dense_diagonal_true) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::SpinUnresolvedONVBasis fock_space {6, 4};
    GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 0};  // 4 alpha 0 beta product SpinUnresolvedONV basis as selected SpinUnresolvedONV basis constructor argument will mimic a spin orbital SpinUnresolvedONV basis
    GQCP::SpinResolvedSelectedONVBasis selected_fock_space {product_fock_space};

    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    // Test the evaluation of the operators with selected SpinUnresolvedONV basis (the reference) versus that of the product SpinUnresolvedONV basis
    auto one_electron_evaluation1 = fock_space.evaluateOperatorDense(h, true);
    auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(h, true);

    auto two_electron_evaluation1 = fock_space.evaluateOperatorDense(g, true);
    auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(g, true);

    auto hamiltonian_evaluation1 = fock_space.evaluateOperatorDense(sq_hamiltonian, true);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(sq_hamiltonian, true);

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


/**
 *  Perform a dense evaluation of a one-, two-electron operator and the Hamiltonian in the SpinUnresolvedONV basis (excluding the diagonal)
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_Dense_diagonal_false) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::SpinUnresolvedONVBasis fock_space {6, 4};
    GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 0};  // 4 alpha 0 beta product SpinUnresolvedONV basis as selected SpinUnresolvedONV basis constructor argument will mimic a spin orbital SpinUnresolvedONV basis
    GQCP::SpinResolvedSelectedONVBasis selected_fock_space {product_fock_space};

    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    // Test the evaluation of the operators with selected SpinUnresolvedONV basis (the reference) versus that of the product SpinUnresolvedONV basis
    auto one_electron_evaluation1 = fock_space.evaluateOperatorDense(h, false);
    auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(h, false);

    auto two_electron_evaluation1 = fock_space.evaluateOperatorDense(g, false);
    auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDense(g, false);

    auto hamiltonian_evaluation1 = fock_space.evaluateOperatorDense(sq_hamiltonian, false);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDense(sq_hamiltonian, false);

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


/**
 *  Perform a sparse evaluation of a one-, two-electron operator and the Hamiltonian in the SpinUnresolvedONV basis (including the diagonal)
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_Sparse_diagonal_true) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::SpinUnresolvedONVBasis fock_space {6, 4};
    GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 0};  // 4 alpha 0 beta product SpinUnresolvedONV basis as selected SpinUnresolvedONV basis constructor argument will mimic a spin orbital SpinUnresolvedONV basis
    GQCP::SpinResolvedSelectedONVBasis selected_fock_space {product_fock_space};

    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    // Test the evaluation of the operators with selected SpinUnresolvedONV basis (the reference) versus that of the product SpinUnresolvedONV basis
    auto one_electron_evaluation1 = GQCP::SquareMatrix<double>(fock_space.evaluateOperatorSparse(h, true));
    auto one_electron_evaluation2 = GQCP::SquareMatrix<double>(selected_fock_space.evaluateOperatorSparse(h, true));

    auto two_electron_evaluation1 = GQCP::SquareMatrix<double>(fock_space.evaluateOperatorSparse(g, true));
    auto two_electron_evaluation2 = GQCP::SquareMatrix<double>(selected_fock_space.evaluateOperatorSparse(g, true));

    auto hamiltonian_evaluation1 = GQCP::SquareMatrix<double>(fock_space.evaluateOperatorSparse(sq_hamiltonian, true));
    auto hamiltonian_evaluation2 = GQCP::SquareMatrix<double>(selected_fock_space.evaluateOperatorSparse(sq_hamiltonian, true));

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


/**
 *  Perform a sparse evaluation of a one-, two-electron operator and the Hamiltonian in the SpinUnresolvedONV basis (excluding the diagonal)
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_Sparse_diagonal_false) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::SpinUnresolvedONVBasis fock_space {6, 4};
    GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 0};  // 4 alpha 0 beta product SpinUnresolvedONV basis as selected SpinUnresolvedONV basis constructor argument will mimic a spin orbital SpinUnresolvedONV basis
    GQCP::SpinResolvedSelectedONVBasis selected_fock_space {product_fock_space};

    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    // Test the evaluation of the operators with selected SpinUnresolvedONV basis (the reference) versus that of the product SpinUnresolvedONV basis
    auto one_electron_evaluation1 = GQCP::SquareMatrix<double>(fock_space.evaluateOperatorSparse(h, false));
    auto one_electron_evaluation2 = GQCP::SquareMatrix<double>(selected_fock_space.evaluateOperatorSparse(h, false));

    auto two_electron_evaluation1 = GQCP::SquareMatrix<double>(fock_space.evaluateOperatorSparse(g, false));
    auto two_electron_evaluation2 = GQCP::SquareMatrix<double>(selected_fock_space.evaluateOperatorSparse(g, false));

    auto hamiltonian_evaluation1 = GQCP::SquareMatrix<double>(fock_space.evaluateOperatorSparse(sq_hamiltonian, false));
    auto hamiltonian_evaluation2 = GQCP::SquareMatrix<double>(selected_fock_space.evaluateOperatorSparse(sq_hamiltonian, false));

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


/**
 *  Evaluate the diagonal of a one-, two-electron operator and the Hamiltonian in the SpinUnresolvedONV basis 
 *  and compare these to the selected CI solutions.
 */
BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_diagonal) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::SpinUnresolvedONVBasis fock_space {6, 4};
    GQCP::SpinResolvedONVBasis product_fock_space {6, 4, 0};  // 4 alpha 0 beta product SpinUnresolvedONV basis as selected SpinUnresolvedONV basis constructor argument will mimic a spin orbital SpinUnresolvedONV basis
    GQCP::SpinResolvedSelectedONVBasis selected_fock_space {product_fock_space};

    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    // Test the evaluation of the operators with selected SpinUnresolvedONV basis (the reference) versus that of the product SpinUnresolvedONV basis
    auto one_electron_evaluation1 = fock_space.evaluateOperatorDiagonal(h);
    auto one_electron_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(h);

    auto two_electron_evaluation1 = fock_space.evaluateOperatorDiagonal(g);
    auto two_electron_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(g);

    auto hamiltonian_evaluation1 = fock_space.evaluateOperatorDiagonal(sq_hamiltonian);
    auto hamiltonian_evaluation2 = selected_fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}


/**
 *  Check the Dense evaluations with diagonal to that of the Dense with the diagonal excluded + the diagonal individually for the Hamiltonian
 */
BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_diagonal_vs_no_diagonal) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::SpinUnresolvedONVBasis fock_space {6, 4};

    GQCP::SquareMatrix<double> hamiltonian = fock_space.evaluateOperatorDense(sq_hamiltonian, true);
    GQCP::SquareMatrix<double> hamiltonian_no_diagonal = fock_space.evaluateOperatorDense(sq_hamiltonian, false);
    GQCP::VectorX<double> hamiltonian_diagonal = fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

    // Test if non-diagonal evaluation and diagonal evaluations are correct
    BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));
}


/**
 *  Perform a matrix-vector product evaluation of a one-, two-electron operator and the Hamiltonian in the SpinUnresolvedONV basis
 *  and compare these to the matrix-vector product of the actual dense evaluations
 */
BOOST_AUTO_TEST_CASE(ONVBasis_EvaluateOperator_MatrixVectorProduct) {

    GQCP::Molecule hchain = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {hchain, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis

    GQCP::SpinUnresolvedONVBasis fock_space {6, 4};

    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    // Generate diagonals for the matvec input
    auto one_electron_diagonal = fock_space.evaluateOperatorDiagonal(h);
    auto two_electron_diagonal = fock_space.evaluateOperatorDiagonal(g);
    auto hamiltonian_diagonal = fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

    // Test the evaluation of the operators with selected SpinUnresolvedONV basis (the reference) versus that of the product SpinUnresolvedONV basis
    auto one_electron_evaluation1 = fock_space.evaluateOperatorMatrixVectorProduct(h, one_electron_diagonal, one_electron_diagonal);
    GQCP::VectorX<double> one_electron_evaluation2 = fock_space.evaluateOperatorDense(h, true) * one_electron_diagonal;

    auto two_electron_evaluation1 = fock_space.evaluateOperatorMatrixVectorProduct(g, two_electron_diagonal, two_electron_diagonal);
    GQCP::VectorX<double> two_electron_evaluation2 = fock_space.evaluateOperatorDense(g, true) * two_electron_diagonal;

    auto hamiltonian_evaluation1 = fock_space.evaluateOperatorMatrixVectorProduct(sq_hamiltonian, hamiltonian_diagonal, hamiltonian_diagonal);
    GQCP::VectorX<double> hamiltonian_evaluation2 = fock_space.evaluateOperatorDense(sq_hamiltonian, true) * hamiltonian_diagonal;

    BOOST_CHECK(one_electron_evaluation1.isApprox(one_electron_evaluation2));
    BOOST_CHECK(two_electron_evaluation1.isApprox(two_electron_evaluation2));
    BOOST_CHECK(hamiltonian_evaluation1.isApprox(hamiltonian_evaluation2));
}

/**
 * Compare the "old" evaluate method with the "new" one that uses the ONVPath API.
 */

BOOST_AUTO_TEST_CASE(ONVBasis_evaluate) {

    // Set up an example molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HChain(5, 0.742, +2);
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto M = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in the orthonormal Löwdin basis

    // Set up the two equivalent ONV bases.
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

    // Check the evaluation of the core Hamiltonian.
    const auto& h_op = sq_hamiltonian.core();

    // Check the dense evaluation.
    const auto h_dense = onv_basis.evaluate<GQCP::SquareMatrix<double>>(h_op, true);  // true: calculate diagonal values
    const auto h_new = onv_basis.evaluate_new<GQCP::SquareMatrix<double>>(h_op, true);

    BOOST_CHECK(h_dense.isApprox(h_new));
}
