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

#define BOOST_TEST_MODULE "OrbitalRotationGenerators_test"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/OrbitalRotationGenerators.hpp"


/**
 *  Check the OrbitalRotationGenerators constructor.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    // Check if the constructor throws upon receiving an argument of wrong dimensions.
    const GQCP::VectorX<double> kappa {4};  // '4' is not a triangular number.
    BOOST_CHECK_THROW(GQCP::OrbitalRotationGenerators generators {kappa}, std::invalid_argument);
}


/**
 *  Check if the `asMatrix` API correctly converts OrbitalRotationGenerators to their matrix representation.
 */
BOOST_AUTO_TEST_CASE(asMatrix) {

    // Set up the raw orbital rotation generators and their matrix representation.
    GQCP::VectorX<double> kappa {3};
    kappa << 1, 2, 3;

    GQCP::SquareMatrix<double> kappa_matrix {3};
    // clang-format off
    kappa_matrix << 0, -1, -2,
                    1,  0, -3,
                    2,  3,  0;
    // clang-format on


    // Check if the wrapper object behaves correctly.
    const GQCP::OrbitalRotationGenerators generators {kappa};
    BOOST_CHECK(generators.asMatrix().isApprox(kappa_matrix));
}


/**
 *  Check if the `FromOccOcc` named constructor behaves correctly.
 */
BOOST_AUTO_TEST_CASE(FromOccOcc) {

    // Initialize a test set of orbital rotation generators for an occupied-occupied orbital space.
    GQCP::VectorX<double> kappa {3};
    kappa << 1, 2, 3;
    const GQCP::OrbitalRotationGenerators occ_occ_generators {kappa};  // 3 (doubly-)occupied spatial orbitals.

    // Construct the total generators in a larger orbital space.
    const auto full_generators = GQCP::OrbitalRotationGenerators::FromOccOcc(occ_occ_generators, 4);  // 4 total spatial orbitals.

    GQCP::SquareMatrix<double> full_kappa_matrix {4};
    // clang-format off
    full_kappa_matrix << 0, -1, -2, 0,
                         1,  0, -3, 0,
                         2,  3,  0, 0,
                         0,  0,  0, 0;
    // clang-format on
    BOOST_CHECK(full_generators.asMatrix().isApprox(full_kappa_matrix));
}
