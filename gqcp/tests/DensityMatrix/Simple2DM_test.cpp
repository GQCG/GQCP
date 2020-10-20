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

#define BOOST_TEST_MODULE "Simple2DM"

#include <boost/test/unit_test.hpp>

#include "DensityMatrix/G2DM.hpp"


/*
 *  MARK: Helper functions
 */

/**
 *  @return a toy 2-DM where
 *      d(i,j,k,l) = l + 2k + 4j + 8i
 */
GQCP::G2DM<double> calculateToy2DMTensor() {

    GQCP::G2DM<double> d {2};

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    auto i_ = static_cast<double>(i);
                    auto j_ = static_cast<double>(j);
                    auto k_ = static_cast<double>(k);
                    auto l_ = static_cast<double>(l);

                    d(i, j, k, l) = l_ + 2 * k_ + 4 * j_ + 8 * i_;
                }
            }
        }
    }

    return d;
};


/*
 *  MARK: Unit tests
 */


/**
 *  Check if the 2-DM `trace` method is correctly implemented, from a manual calculation.
 */
BOOST_AUTO_TEST_CASE(trace) {

    const auto d = calculateToy2DMTensor();

    BOOST_CHECK(std::abs(d.trace() - 30.0) < 1.0e-12);
}


/**
 *  Check if the 2-DM `reduce` method is correctly implemented, from a manual calculation.
 */
BOOST_AUTO_TEST_CASE(reduce) {

    const auto d = calculateToy2DMTensor();

    // Set up the reference result.
    GQCP::G1DM<double> D_ref = GQCP::G1DM<double>::Zero(2);

    // clang-format off
    D_ref <<  3, 11,
             19, 27;
    // clang-format on

    BOOST_CHECK(D_ref.isApprox(d.reduce(), 1.0e-12));
}
