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
#include "QCModel/CI/LinearExpansion.hpp"


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
 *  Test if the SpinUnresolvedONV basis correctly counts the number of coupling ONVs with larger address for a given SpinUnresolvedONV.
 */
BOOST_AUTO_TEST_CASE(coupling_count) {

    const GQCP::SpinUnresolvedONVBasis onv_basis {5, 3};
    GQCP::SpinUnresolvedONV onv = onv_basis.constructONVFromAddress(3);  // 01110

    // We only count couplings with larger addresses.
    BOOST_CHECK(onv_basis.countOneElectronCouplings(onv) == 3);      // 11100, 11010, 10110
    BOOST_CHECK(onv_basis.countTwoElectronCouplings(onv) == 3 + 3);  // 11100, 11010, 10110, 11001, 10101, 10011

    onv = onv_basis.constructONVFromAddress(0);  // 00111

    BOOST_CHECK(onv_basis.countOneElectronCouplings(onv) == 6);
    BOOST_CHECK(onv_basis.countTwoElectronCouplings(onv) == 6 + 3);  // All of them.


    // Test whether the total count matches that of individual counts of all ONVs in the SpinUnresolvedONV basis.
    const GQCP::SpinUnresolvedONVBasis onv_basis2 {16, 8};

    size_t coupling_count1 = 0;
    size_t coupling_count2 = 0;
    onv = onv_basis2.constructONVFromAddress(0);  // ONV with address 0.
    for (size_t I = 0; I < onv_basis2.dimension(); I++) {
        if (I > 0) {
            onv_basis2.transformONVToNextPermutation(onv);
        }
        coupling_count1 += onv_basis2.countOneElectronCouplings(onv);
        coupling_count2 += onv_basis2.countTwoElectronCouplings(onv);
    }

    BOOST_CHECK(2 * coupling_count1 == onv_basis2.countTotalOneElectronCouplings());
    BOOST_CHECK(2 * coupling_count2 == onv_basis2.countTotalTwoElectronCouplings());
}


/**
 *  In this test we iterate over the entire SpinUnresolvedONV basis using `transformONVToNextPermutation`. We test if the address is correct using `address`.
 */
BOOST_AUTO_TEST_CASE(ONV_address_transform_fullspace) {

    // Here we will test a full permutation through a SpinUnresolvedONV basis of K = 15, N = 5.
    const GQCP::SpinUnresolvedONVBasis onv_basis {15, 5};

    // Retrieve the first SpinUnresolvedONV of the SpinUnresolvedONV basis.
    GQCP::SpinUnresolvedONV onv_test = onv_basis.constructONVFromAddress(0);

    const size_t dimension_fock_space = 3003;
    bool is_correct = true;  // A variable that is updated to false if an unexpected result occurs.

    // Iterate through the SpinUnresolvedONV basis in reverse lexicographical order and test whether address matches.
    for (size_t i = 0; i < dimension_fock_space; i++) {

        // Test the address.
        if (i != onv_basis.addressOf(onv_test)) {
            is_correct = false;
        }

        // Transform the given SpinUnresolvedONV to the next SpinUnresolvedONV.
        if (i < dimension_fock_space - 1) {
            onv_basis.transformONVToNextPermutation(onv_test);
        }
    }

    // Checks if no unexpected results occured in a full iteration
    BOOST_CHECK(is_correct);
}


/**
 *  Test whether the SpinUnresolvedONV basis attributes the correct address to a given SpinUnresolvedONV.
 */
BOOST_AUTO_TEST_CASE(ONVBasis_addressOf) {

    const GQCP::SpinUnresolvedONVBasis onv_basis {6, 3};

    // The address of the string "010011" (19) should be 4.
    const GQCP::SpinUnresolvedONV onv {6, 3, 19};

    BOOST_CHECK_EQUAL(onv_basis.addressOf(onv), 4);
}


/**
 *  Test the `transformONVToNextPermutation` API for manually chosen ONVs.
 */
BOOST_AUTO_TEST_CASE(transformONVToNextPermutation) {

    const GQCP::SpinUnresolvedONVBasis onv_basis {5, 3};
    // K = 5, N = 3 <-> "00111"
    GQCP::SpinUnresolvedONV onv = onv_basis.constructONVFromAddress(0);
    // The lexical permutations are: "00111" (7), "01011" (11), "01101" (13), "01110" (14), etc.

    // Check the permutations one after the other.

    onv_basis.transformONVToNextPermutation(onv);  // "01011" (11)
    BOOST_CHECK_EQUAL(onv.unsignedRepresentation(), 11);
    const std::vector<size_t> ref_indices1 {0, 1, 3};
    BOOST_CHECK(ref_indices1 == onv.occupiedIndices());

    onv_basis.transformONVToNextPermutation(onv);  // "01101" (13)
    BOOST_CHECK_EQUAL(onv.unsignedRepresentation(), 13);
    const std::vector<size_t> ref_indices2 {0, 2, 3};
    BOOST_CHECK(ref_indices2 == onv.occupiedIndices());

    onv_basis.transformONVToNextPermutation(onv);  // "01110" (14)
    BOOST_CHECK_EQUAL(onv.unsignedRepresentation(), 14);
    const std::vector<size_t> ref_indices3 {1, 2, 3};
    BOOST_CHECK(ref_indices3 == onv.occupiedIndices());
}


/**
 *  Check if the dense evaluation of a one-electron operator in a spin-unresolved ONV basis matches the evaluation in a spin-resolved selected ONV basis with only alpha electrons.
 * 
 *  The test system is a H6(2+)-chain with internuclear separation of 0.742 (a.u.) in an STO-3G basis.
 */
BOOST_AUTO_TEST_CASE(evaluate_one_electron_operator_dense) {

    // Set up an example molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(6, 0.742, +2);
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto M = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);

    // Set up the two equivalent ONV bases.
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};
    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {GQCP::SpinResolvedONVBasis(M, N, 0)};  // No beta electrons, to mimic a spin-unresolved case.


    // Check the evaluation of the core Hamiltonian. We'll have to convert the restricted operator to a generalized operator in order to use the semantically correct APIs.
    const auto& h = sq_hamiltonian.core();
    const auto h_alpha_generalized = GQCP::ScalarGSQOneElectronOperator<double>::FromUnrestrictedComponent(h.alpha());

    // Check the dense evaluation.
    const auto h_dense = onv_basis.evaluateOperatorDense(h_alpha_generalized);
    const auto h_dense_selected = selected_onv_basis.evaluateOperatorDense(h);
    BOOST_CHECK(h_dense.isApprox(h_dense_selected, 1.0e-12));
}


/**
 *  Check if the dense evaluation of a two-electron operator in a spin-unresolved ONV basis matches the evaluation in a spin-resolved selected ONV basis with only alpha electrons.
 * 
 *  The test system is a H6(2+)-chain with internuclear separation of 0.742 (a.u.) in an STO-3G basis.
 */
BOOST_AUTO_TEST_CASE(evaluate_two_electron_operator_dense) {

    // Set up an example molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(6, 0.742, +2);
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto M = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);

    // Set up the two equivalent ONV bases.
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};
    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {GQCP::SpinResolvedONVBasis(M, N, 0)};  // No beta electrons, to mimic a spin-unresolved case.


    // Check the evaluation of the two-electron part of the Hamiltonian. We'll have to convert the restricted operator to a generalized operator in order to use the semantically correct APIs.
    const auto& g = sq_hamiltonian.twoElectron();
    const auto g_alpha_generalized = GQCP::ScalarGSQTwoElectronOperator<double>::FromUnrestrictedComponent(g.alphaAlpha());

    // Check the dense evaluation.
    const auto g_dense = onv_basis.evaluateOperatorDense(g_alpha_generalized);
    const auto g_dense_selected = selected_onv_basis.evaluateOperatorDense(g);
    BOOST_CHECK(g_dense.isApprox(g_dense_selected, 1.0e-12));
}


/**
 *  Check if the diagonal of the matrix representation of a generalized Hamiltonian is equal to the diagonal that is calculated through a specialized routine.
 * 
 *  The test system is a H6(2+)-chain with internuclear separation of 0.742 (a.u.) in an STO-3G basis.
 */
BOOST_AUTO_TEST_CASE(generalized_hamiltonian_diagonal) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    const auto hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(spinor_basis, molecule);
    const auto K = hamiltonian.numberOfOrbitals();

    // Set up the full spin-unresolved ONV basis.
    const GQCP::SpinUnresolvedONVBasis onv_basis {K, molecule.numberOfElectrons()};

    // Determine the Hamiltonian matrix and the diagonal through a specialized routine, and check if they match.
    const auto dense_matrix = onv_basis.evaluateOperatorDense(hamiltonian);
    const auto diagonal_specialized = onv_basis.evaluateOperatorDiagonal(hamiltonian);

    BOOST_CHECK(diagonal_specialized.isApprox(dense_matrix.diagonal(), 1.0e-12));
}


/**
 *  Check if the sparse evaluation of a one-electron operator in a spin-unresolved ONV basis matches the evaluation in a spin-resolved selected ONV basis with only alpha electrons.
 * 
 *  The test system is a H6(2+)-chain with internuclear separation of 0.742 (a.u.) in an STO-3G basis.
 */
BOOST_AUTO_TEST_CASE(evaluate_one_electron_operator_sparse) {

    // Set up an example molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HChain(6, 0.742, +2);
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto M = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);

    // Set up the two equivalent ONV bases.
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};
    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {GQCP::SpinResolvedONVBasis(M, N, 0)};  // No beta electrons, to mimic a spin-unresolved case.


    // Check the evaluation of the core Hamiltonian. We'll have to convert the restricted operator to a generalized operator in order to use the semantically correct APIs.
    const auto& h = sq_hamiltonian.core();
    const auto h_alpha_generalized = GQCP::ScalarGSQOneElectronOperator<double>::FromUnrestrictedComponent(h.alpha());

    // Check the sparse evaluation.
    const auto h_sparse = onv_basis.evaluateOperatorSparse(h_alpha_generalized);
    const auto h_sparse_selected = selected_onv_basis.evaluateOperatorSparse(h);
    BOOST_CHECK(h_sparse.isApprox(h_sparse_selected, 1.0e-12));
}


/**
 *  Check if the sparse evaluation of a two-electron operator in a spin-unresolved ONV basis matches the evaluation in a spin-resolved selected ONV basis with only alpha electrons.
 * 
 *  The test system is a H6(2+)-chain with internuclear separation of 0.742 (a.u.) in an STO-3G basis.
 */
BOOST_AUTO_TEST_CASE(evaluate_two_electron_operator_sparse) {

    // Set up an example molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HChain(6, 0.742, +2);
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto M = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In the orthonormal Löwdin basis.

    // Set up the two equivalent ONV bases.
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};
    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {GQCP::SpinResolvedONVBasis(M, N, 0)};  // No beta electrons, to mimic a spin-unresolved case.


    // Check the evaluation of the two-electron part of the Hamiltonian. We'll have to convert the restricted operator to a generalized operator in order to use the semantically correct APIs.
    const auto& g = sq_hamiltonian.twoElectron();
    const auto g_alpha_generalized = GQCP::ScalarGSQTwoElectronOperator<double>::FromUnrestrictedComponent(g.alphaAlpha());

    // Check the sparse evaluation.
    const auto g_sparse = onv_basis.evaluateOperatorDense(g_alpha_generalized);
    const auto g_sparse_selected = selected_onv_basis.evaluateOperatorDense(g);
    BOOST_CHECK(g_sparse.isApprox(g_sparse_selected, 1.0e-12));
}


/**
 *  Check if the matrix-vector product through a direct evaluation (i.e. through the dense Hamiltonian matrix representation) and the specialized implementation are equal.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a spin-unresolved FCI dimension of 1001.
 */
BOOST_AUTO_TEST_CASE(generalized_dense_vs_matvec) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    const auto hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(spinor_basis, molecule);
    const auto M = hamiltonian.numberOfOrbitals();

    // Set up the full spin-unresolved ONV basis.
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, molecule.numberOfElectrons()};

    // Determine the Hamiltonian matrix and let it act on a random linear expansion.
    const auto linear_expansion = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::Random(onv_basis);
    const auto H_dense = onv_basis.evaluateOperatorDense(hamiltonian);
    const GQCP::VectorX<double> direct_mvp = H_dense * linear_expansion.coefficients();  // mvp: matrix-vector-product

    // Determine the specialized matrix-vector product and check if they are equal.
    const auto specialized_mvp = onv_basis.evaluateOperatorMatrixVectorProduct(hamiltonian, linear_expansion.coefficients());

    BOOST_CHECK(specialized_mvp.isApprox(direct_mvp, 1.0e-08));
}


/*
 *  MARK: Tests for legacy code
 */

/**
 *  Test if the shift in address and orbital indices correspond to the correct solutions
 *  
 *  This method does not alter the ONVs in any way.
 */
BOOST_AUTO_TEST_CASE(iterateToNextUnoccupiedOrbital) {

    GQCP::SpinUnresolvedONVBasis onv_basis {5, 3};
    GQCP::SpinUnresolvedONV onv = onv_basis.constructONVFromAddress(3);  // 01110

    size_t address_shift = 0;
    // test shift if we annihilate one electron and start from orbital index 2
    size_t e = 1;  // count starts at 1 (translates to orbital index 2)
    size_t q = 2;  // index starts at orbital index 2

    //  In this instance electron weights at index 2 and 3 should be shifted.
    //  Initial weight contributions were 1 and 1 respectively,
    //  these should be shifted to 2 and 3 respectively, the difference is 1 and 2 respectively.
    //  The total shift is thus 3
    onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address_shift, q, e);

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
    onv_basis.shiftUntilNextUnoccupiedOrbital<2>(onv, address_shift, q, e);

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

    GQCP::SpinUnresolvedONVBasis onv_basis {5, 3};
    GQCP::SpinUnresolvedONV onv = onv_basis.constructONVFromAddress(3);  // 01110

    size_t address_shift = 0;
    int sign = 1;
    // test shift if we annihilate one electron and start from orbital index 2
    size_t e = 1;  // count starts at 1 (translates to orbital index 2)
    size_t q = 2;  // index starts at orbital index 2

    //  In this instance electron weights at index 2 and 3 should be shifted.
    //  Initial weight contributions were 1 and 1 respectively,
    //  these should be shifted to 2 and 3 respectively, the difference is 1 and 2 respectively.
    //  The total shift is thus 3, and the sign should remain the same (flips twice).
    onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address_shift, q, e, sign);

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
    onv_basis.shiftUntilNextUnoccupiedOrbital<2>(onv, address_shift, q, e, sign);

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

    GQCP::SpinUnresolvedONVBasis onv_basis {5, 3};
    GQCP::SpinUnresolvedONV onv = onv_basis.constructONVFromAddress(6);  // 10110

    size_t address_shift = 0;
    int sign = 1;

    // test shift if we plan on creating one electron and start from orbital index 2
    size_t e = 0;  // count starts at 0 (translates to orbital index 1)
    size_t q = 1;  // index starts at orbital index 1

    //  Index 1 is occupied an thus its weight shall shift if we create before an electron on a smaller index.
    //  In this instance the electron weight at index 1 should be shifted.
    //  Initial weight contribution is 1,
    //  This should be shifted to 0, the difference is 1.
    //  The sign changes (flips once).
    onv_basis.shiftUntilPreviousUnoccupiedOrbital<1>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == -1);
    BOOST_CHECK(e == -1);
    BOOST_CHECK(q == 0);
    BOOST_CHECK(sign == -1);

    sign = 1;
    address_shift = 0;
    onv = onv_basis.constructONVFromAddress(9);  // 11100
    // Test shift if we plan on creating two electrons and start from orbital index 2.
    e = 0;  // count starts at 1 (translates to orbital index 2)
    q = 2;  // index starts at orbital index 2

    //  Index 2 is occupied an thus its weight shall shift if we create before an electron on a smaller index.
    //  In this instance the electron weight at index 2 should be shifted.
    //  Initial weight contribution is 2,
    //  This should be shifted to 0, the difference is 2.
    //  The sign changes (flips once).
    onv_basis.shiftUntilPreviousUnoccupiedOrbital<2>(onv, address_shift, q, e, sign);

    BOOST_CHECK(address_shift == -2);
    BOOST_CHECK(e == -1);
    BOOST_CHECK(q == 1);
    BOOST_CHECK(sign == -1);
}
