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

#define BOOST_TEST_MODULE "OperatorString"

#include <boost/test/unit_test.hpp>

#include "ONVBasis/SpinResolvedOperatorString.hpp"
#include "ONVBasis/SpinUnresolvedOperatorString.hpp"
#include "QuantumChemical/Spin.hpp"


/**
 *  Check whether the SpinUnresolvedOperatorString constructor works when a `SpinUnresolvedONV` is given as parameter.
 */
BOOST_AUTO_TEST_CASE(test_unresolved_onv_constructor) {

    const GQCP::SpinUnresolvedONV onv {4, 2, 5};  // "1010"
    const GQCP::SpinUnresolvedOperatorString operator_string = GQCP::SpinUnresolvedOperatorString::FromONV(onv);

    const std::vector<size_t> correct_indices = {0, 2};

    BOOST_CHECK_EQUAL_COLLECTIONS(operator_string.operatorIndices().begin(), operator_string.operatorIndices().end(), correct_indices.begin(), correct_indices.end());
}


/**
 *  Check whether the SpinResolvedOperatorString constructor works.
 */
BOOST_AUTO_TEST_CASE(test_resolved_constructor) {

    // Initialize the vector containing the indices we want to use for our operator string.
    std::vector<size_t> indices = {1, 4, 7, 8};

    // Initialize the vector containing the indices we want to use for our operator string.
    // However, the dimension of this index vector doesn't match the spin vector dimension and should result in an exception being thrown.
    std::vector<size_t> indices_wrong = {4, 7, 8};

    // Initialize the vector containing the spin components we want to use for our operator string indices.
    std::vector<GQCP::Spin> spins = {GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha, GQCP::Spin::beta};

    // Initialize the vector containing the spin components we want to use for our operator string indices.
    // However, the dimension of this spin vector doesn't match the index vector dimension and should result in an exception being thrown.
    std::vector<GQCP::Spin> spins_wrong = {GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha};

    // Check the constructor.
    BOOST_CHECK_THROW(GQCP::SpinResolvedOperatorString operator_string_wrong(indices_wrong, spins), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::SpinResolvedOperatorString operator_string_wrong_2(indices, spins_wrong), std::invalid_argument);
    BOOST_CHECK_NO_THROW(GQCP::SpinResolvedOperatorString operator_string_correct(indices, spins));
}


/**
 * 
 * 
 */
BOOST_AUTO_TEST_CASE(test_resolved_onv_constructor) {

    const GQCP::SpinUnresolvedONV onv_alpha {4, 2, 5};  //  "1010"
    const GQCP::SpinUnresolvedONV onv_beta {4, 2, 10};  //  "0101"
    const GQCP::SpinResolvedONV onv {onv_alpha, onv_beta};

    std::cout << "\nonv_alpha indices" << std::endl;
    for (const size_t& i : onv.onv(GQCP::Spin::alpha).occupiedIndices()) {
        std::cout << i << "\t";
    }
    std::cout << "\nonv_beta indices" << std::endl;
    for (const size_t& i : onv.onv(GQCP::Spin::beta).occupiedIndices()) {
        std::cout << i << "\t";
    }

    const GQCP::SpinResolvedOperatorString operator_string = GQCP::SpinResolvedOperatorString::FromONV(onv);
    const auto operator_indices = operator_string.operatorIndices();
    const auto operator_spins = operator_string.operatorSpins();
    const std::vector<size_t> correct_indices = {0, 2, 1, 3};                                                                  // First the indices of the alpha operator string, then those of the beta operator string.
    const std::vector<GQCP::Spin> correct_spins = {GQCP::Spin::alpha, GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::beta};  // First the alpha spins, then the beta spins.

    BOOST_CHECK_EQUAL_COLLECTIONS(operator_indices.begin(), operator_indices.end(), correct_indices.begin(), correct_indices.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(operator_spins.begin(), operator_spins.end(), correct_spins.begin(), correct_spins.end());
}


/**
 *  Check whether the SpinUnresolvedOperatorString stores its information correctly.
 */
BOOST_AUTO_TEST_CASE(test_unresolved_indices) {

    // Initialize the vector containing the indices we want to use for our operator string.
    std::vector<size_t> indices = {1, 4, 7, 8};

    // Initialize the SpinUnresolvedOperatorString.
    const GQCP::SpinUnresolvedOperatorString operator_string {indices};

    // Get the indices from the operator string.
    const auto operator_string_indices = operator_string.operatorIndices();

    // Check whether the indices are stored correctly.
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), operator_string_indices.begin(), operator_string_indices.end());
}


/**
 *  Check whether the SpinUnresolvedOperatorString`s `isZero()` check correctly checks whether the operator string will inherently equal zero when applied to any ONV.
 */
BOOST_AUTO_TEST_CASE(test_unresolved_isZero) {

    // Initialize the vector containing the indices we want to use for our operator string.
    // One vector does contain duplicate indices, one does not.
    std::vector<size_t> indices_1 = {1, 4, 7, 8};
    std::vector<size_t> indices_2 = {1, 4, 7, 1};


    // Initialize the SpinUnresolvedOperatorStrings corresponding to the indices.
    const GQCP::SpinUnresolvedOperatorString operator_string_1 {indices_1};
    const GQCP::SpinUnresolvedOperatorString operator_string_2 {indices_2};

    // Compare the result of the `isZero()` check with the expected boolean value (false for the operator string without duplicate indices, true for the one with duplicate indices).
    BOOST_CHECK_EQUAL(operator_string_1.isZero(), false);
    BOOST_CHECK_EQUAL(operator_string_2.isZero(), true);
}


/**
 *  Check whether the SpinResolvedOperatorString stores its information correctly.
 */
BOOST_AUTO_TEST_CASE(test_resolved_indices_spins) {

    // Initialize the vector containing the indices we want to use for our operator string.
    std::vector<size_t> indices = {1, 4, 7, 8};

    // Initialize the vector containing the spin components we want to use for our operator string indices.
    std::vector<GQCP::Spin> spins = {GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha, GQCP::Spin::beta};

    // Initialize the SpinUnresolvedOperatorString.
    const GQCP::SpinResolvedOperatorString operator_string {indices, spins};

    // Get the indices and spins from the operator string.
    const auto operator_string_indices = operator_string.operatorIndices();
    const auto operator_string_spins = operator_string.operatorSpins();

    // Check whether the indices and spins are stored correctly.
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), operator_string_indices.begin(), operator_string_indices.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(spins.begin(), spins.end(), operator_string_spins.begin(), operator_string_spins.end());
}


/**
 *  Check whether the SpinResolvedOperatorString can be correctly `SpinResolved` in its alpha and beta component, with the correct phase factors.
 */
BOOST_AUTO_TEST_CASE(test_spin_resolve) {

    // Initialize the vector containing the indices we want to use for our operator string.
    std::vector<size_t> indices = {0, 1, 2, 3, 4};

    // Initialize the vector containing the spin components we want to use for our operator string indices.
    std::vector<GQCP::Spin> spins = {GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha};

    // Initialize the SpinUnresolvedOperatorString.
    const auto operator_string = GQCP::SpinResolvedOperatorString(indices, spins);

    // `Spin resolve` the operator string, meaning we split it in it's alpha and beta components, with the appropriate phase factors.
    const auto spin_resolved_operator_string = operator_string.spinResolve();

    // Get the indices from each new operator string.
    const auto alpha_operator_string_indices = spin_resolved_operator_string.alpha().operatorIndices();
    const auto beta_operator_string_indices = spin_resolved_operator_string.beta().operatorIndices();

    // Define references to check the results.
    // The phase factor of the alpha string should be -1:
    //
    //      a_a a_b a_a a_b a_a (initial string) & p = 1 (initial phase factor)
    //      => a_a a_a a_b a_b a_a (first swap) & p = -1 (phase factor after swap 1)
    //      => a_a a_a a_b a_a a_b (second swap) & p = 1 (phase factor after swap 2)
    //      => a_a a_a a_a a_b a_b (third swap) & p = -1 (phase factor after swap 3)
    //
    // Now we have:
    //
    //      (-1) a_a a_a a_a a_b a_b |vac_a> |vac_b>

    const auto reference_phase_factor = -1;

    // The alpha and beta indices respectively can easily be seen from the initial pairs.
    std::vector<size_t> alpha_indices = {0, 2, 4};
    std::vector<size_t> beta_indices = {1, 3};

    // Test the calculated values against the reference values.
    BOOST_CHECK_EQUAL_COLLECTIONS(alpha_indices.begin(), alpha_indices.end(), alpha_operator_string_indices.begin(), alpha_operator_string_indices.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(beta_indices.begin(), beta_indices.end(), beta_operator_string_indices.begin(), beta_operator_string_indices.end());
    BOOST_CHECK_EQUAL(reference_phase_factor, spin_resolved_operator_string.alpha().phaseFactor());
}


/**
 *  Check whether the SpinResolvedOperatorString can be correctly `SpinResolved` in its alpha and beta component, with the correct phase factors.
 */
BOOST_AUTO_TEST_CASE(test_spin_resolve_2) {

    // Initialize the vector containing the indices we want to use for our operator string.
    std::vector<size_t> indices = {0, 1, 2, 3, 4, 5};

    // Initialize the vector containing the spin components we want to use for our operator string indices.
    std::vector<GQCP::Spin> spins = {GQCP::Spin::beta, GQCP::Spin::beta, GQCP::Spin::beta, GQCP::Spin::alpha, GQCP::Spin::alpha, GQCP::Spin::alpha};

    // Initialize the SpinUnresolvedOperatorString.
    const auto operator_string = GQCP::SpinResolvedOperatorString(indices, spins);

    // `Spin resolve` the operator string, meaning we split it in it's alpha and beta components, with the appropriate phase factors.
    const auto spin_resolved_operator_string = operator_string.spinResolve();

    // Get the indices from each new operator string.
    const auto alpha_operator_string_indices = spin_resolved_operator_string.alpha().operatorIndices();
    const auto beta_operator_string_indices = spin_resolved_operator_string.beta().operatorIndices();

    // Define references to check the results.
    // The phase factor of the alpha string should be -1:
    //
    //      a_b a_b a_b a_a a_a a_a (initial string) & p = 1 (initial phase factor)
    //      => a_3a has three betas in front of it, p = (-1)^3
    //      => a_4a has three betas in front of it, p = (-1)^3
    //      => a_5a has three betas in front of it, p = (-1)^3
    //
    // As a total phase factor for the alpha string we now have p = (-1)^9 = -1

    const auto reference_phase_factor = -1;

    // The alpha and beta indices respectively can easily be seen from the initial pairs.
    std::vector<size_t> alpha_indices = {3, 4, 5};
    std::vector<size_t> beta_indices = {0, 1, 2};

    // Test the calculated values against the reference values.
    BOOST_CHECK_EQUAL_COLLECTIONS(alpha_indices.begin(), alpha_indices.end(), alpha_operator_string_indices.begin(), alpha_operator_string_indices.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(beta_indices.begin(), beta_indices.end(), beta_operator_string_indices.begin(), beta_operator_string_indices.end());
    BOOST_CHECK_EQUAL(reference_phase_factor, spin_resolved_operator_string.alpha().phaseFactor());
}


/**
 *  Check whether the SpinResolvedOperatorString can be correctly `SpinResolved` in its alpha and beta component, with the correct phase factors.
 */
BOOST_AUTO_TEST_CASE(test_spin_resolve_3) {

    // Initialize the vector containing the indices we want to use for our operator string.
    std::vector<size_t> indices = {7, 8, 4, 1, 2, 4, 3, 7};

    // Initialize the vector containing the spin components we want to use for our operator string indices.
    std::vector<GQCP::Spin> spins = {GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha, GQCP::Spin::beta};

    // Initialize the SpinUnresolvedOperatorString.
    const auto operator_string = GQCP::SpinResolvedOperatorString(indices, spins);

    // `Spin resolve` the operator string, meaning we split it in it's alpha and beta components, with the appropriate phase factors.
    const auto spin_resolved_operator_string = operator_string.spinResolve();

    // Get the indices from each new operator string.
    const auto alpha_operator_string_indices = spin_resolved_operator_string.alpha().operatorIndices();
    const auto beta_operator_string_indices = spin_resolved_operator_string.beta().operatorIndices();

    // Define references to check the results.
    // The phase factor of the alpha string should be 1:
    //
    //      a_a a_b a_a a_b a_a a_b a_a a_b (initial string) & p = 1 (initial phase factor)
    //      => a_7a has zero betas in front of it, p = 1
    //      => a_4a has one beta in front of it, p = -1
    //      => a_2a has two betas in front of it, p = (-1)^2
    //      => a_3a has three betas in front of it, p = (-1)^3
    //
    // As a total phase factor for the alpha string we now have p = (-1)^6 = 1

    const auto reference_phase_factor = 1;

    // The alpha and beta indices respectively can easily be seen from the initial pairs.
    std::vector<size_t> alpha_indices = {7, 4, 2, 3};
    std::vector<size_t> beta_indices = {8, 1, 4, 7};

    // Test the calculated values against the reference values.
    BOOST_CHECK_EQUAL_COLLECTIONS(alpha_indices.begin(), alpha_indices.end(), alpha_operator_string_indices.begin(), alpha_operator_string_indices.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(beta_indices.begin(), beta_indices.end(), beta_operator_string_indices.begin(), beta_operator_string_indices.end());
    BOOST_CHECK_EQUAL(reference_phase_factor, spin_resolved_operator_string.alpha().phaseFactor());
}
