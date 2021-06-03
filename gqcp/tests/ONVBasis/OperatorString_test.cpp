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
 *  Check whether the SpinUnresolvedOperatorString constructor works.
 */
BOOST_AUTO_TEST_CASE(test_unresolved_constructor) {

    // Initialize the vector which we want to use for our operator string.
    std::vector<size_t> indices = {1, 4, 7, 8};

    // Check the constructor.
    BOOST_CHECK_NO_THROW(GQCP::SpinUnresolvedOperatorString operator_string {indices});
}

/**
 *  Check whether the SpinResolvedOperatorString constructor works.
 */
BOOST_AUTO_TEST_CASE(test_resolved_constructor) {

    // Initialize the vector which we want to use for our operator string indices.
    std::vector<size_t> indices = {1, 4, 7, 8};

    // Initialize the vector which we want to use for our operator string indices.
    std::vector<size_t> indices_wrong = {4, 7, 8};

    // Initialize the vector which we want to use for our operator string indices.
    std::vector<GQCP::Spin> spins = {GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha, GQCP::Spin::beta};

    // Check the constructor.
    BOOST_CHECK_THROW(GQCP::SpinResolvedOperatorString operator_string_wrong(indices_wrong, spins), std::invalid_argument);
    BOOST_CHECK_NO_THROW(GQCP::SpinResolvedOperatorString operator_string_correct(indices, spins));
}

/**
 *  Check whether the SpinUnresolvedOperatorString stores its information correctly.
 */
BOOST_AUTO_TEST_CASE(test_unresolved_indices) {

    // Initialize the vector which we want to use for our operator string.
    std::vector<size_t> indices = {1, 4, 7, 8};

    // initialize the SpinUnresolvedOperatorString.
    const auto operator_string = GQCP::SpinUnresolvedOperatorString {indices};

    // Get the indices from the operator string.
    const auto operator_string_indices = operator_string.operatorIndices();

    // Check whether the indices are stored correctly.
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), operator_string_indices.begin(), operator_string_indices.end());
}

/**
 *  Check whether the SpinResolvedOperatorString stores its information correctly.
 */
BOOST_AUTO_TEST_CASE(test_resolved_indices_spins) {

    // Initialize the vector which we want to use for our operator string.
    std::vector<size_t> indices = {1, 4, 7, 8};

    // Initialize the vector which we want to use for our operator string indices.
    std::vector<GQCP::Spin> spins = {GQCP::Spin::alpha, GQCP::Spin::beta, GQCP::Spin::alpha, GQCP::Spin::beta};

    // initialize the SpinUnresolvedOperatorString.
    const auto operator_string = GQCP::SpinResolvedOperatorString(indices, spins);

    // Get the indices and spins from the operator string.
    const auto operator_string_indices = operator_string.operatorIndices();
    const auto operator_string_spins = operator_string.operatorSpins();

    // Check whether the indices are stored correctly.
    BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), operator_string_indices.begin(), operator_string_indices.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(spins.begin(), spins.end(), operator_string_spins.begin(), operator_string_spins.end());
}
