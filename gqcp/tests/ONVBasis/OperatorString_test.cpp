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

#define BOOST_TEST_MODULE "ONVPath"

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
