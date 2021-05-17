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

#include "Basis/Transformations/GOrbitalRotationGenerators.hpp"
#include "Basis/Transformations/ROrbitalRotationGenerators.hpp"
#include "Basis/Transformations/UOrbitalRotationGenerators.hpp"
#include "Basis/Transformations/UOrbitalRotationGeneratorsComponent.hpp"


/**
 *  Check the ROrbitalRotationGenerators constructor.
 */
BOOST_AUTO_TEST_CASE(restrictedConstructor) {

    // Check if the constructor throws upon receiving an argument of wrong dimensions.
    const GQCP::VectorX<double> kappa {4};  // '4' is not a triangular number.
    BOOST_CHECK_THROW(GQCP::ROrbitalRotationGenerators<double> generators {kappa}, std::invalid_argument);
}

/**
 *  Check the UOrbitalRotationGeneratorsComponent constructor.
 */
BOOST_AUTO_TEST_CASE(unrestrictedComponentConstructor) {

    // Check if the constructor throws upon receiving an argument of wrong dimensions.
    const GQCP::VectorX<double> kappa {4};  // '4' is not a triangular number.
    BOOST_CHECK_THROW(GQCP::UOrbitalRotationGeneratorsComponent<double> generators {kappa}, std::invalid_argument);
}

/**
 *  Check the GOrbitalRotationGenerators constructor.
 */
BOOST_AUTO_TEST_CASE(generalizedConstructor) {

    // Check if the constructor throws upon receiving an argument of wrong dimensions.
    const GQCP::VectorX<double> kappa {4};  // '4' is not a triangular number.
    BOOST_CHECK_THROW(GQCP::GOrbitalRotationGenerators<double> generators {kappa}, std::invalid_argument);
}


/**
 *  Check if the `asMatrix` API correctly converts ROrbitalRotationGenerators to their matrix representation.
 */
BOOST_AUTO_TEST_CASE(restrictedAsMatrix) {

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
    const GQCP::ROrbitalRotationGenerators<double> generators {kappa};
    BOOST_CHECK(generators.asMatrix().isApprox(kappa_matrix));
}


/**
 *  Check if the `asMatrix` API correctly converts UOrbitalRotationGenerators to their matrix representation.
 */
BOOST_AUTO_TEST_CASE(unrestrictedAsMatrix) {

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
    const GQCP::UOrbitalRotationGeneratorsComponent<double> alpha_generators {kappa};
    const GQCP::UOrbitalRotationGeneratorsComponent<double> beta_generators {kappa};

    const GQCP::UOrbitalRotationGenerators<double> generators {alpha_generators, beta_generators};

    BOOST_CHECK(generators.asMatrix().alpha().isApprox(kappa_matrix));
    BOOST_CHECK(generators.asMatrix().beta().isApprox(kappa_matrix));
}


/**
 *  Check if the `asMatrix` API correctly converts GOrbitalRotationGenerators to their matrix representation.
 */
BOOST_AUTO_TEST_CASE(generalizedAsMatrix) {

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
    const GQCP::GOrbitalRotationGenerators<double> generators {kappa};
    BOOST_CHECK(generators.asMatrix().isApprox(kappa_matrix));
}


/**
 *  Check if the `FromOccOcc` named constructor behaves correctly for the restricted case.
 */
BOOST_AUTO_TEST_CASE(restrictedFromOccOcc) {

    // Initialize a test set of orbital rotation generators for an occupied-occupied orbital space.
    GQCP::VectorX<double> kappa {3};
    kappa << 1, 2, 3;
    const GQCP::ROrbitalRotationGenerators<double> occ_occ_generators {kappa};  // 3 (doubly-)occupied spatial orbitals.

    // Construct the total generators in a larger orbital space.
    const auto full_generators = GQCP::ROrbitalRotationGenerators<double>::FromOccOcc(occ_occ_generators, 4);  // 4 total spatial orbitals.

    GQCP::SquareMatrix<double> full_kappa_matrix {4};
    // clang-format off
    full_kappa_matrix << 0, -1, -2, 0,
                         1,  0, -3, 0,
                         2,  3,  0, 0,
                         0,  0,  0, 0;
    // clang-format on
    BOOST_CHECK(full_generators.asMatrix().isApprox(full_kappa_matrix));
}


/**
 *  Check if the `FromOccOcc` named constructor behaves correctly for the unrestricted case.
 */
BOOST_AUTO_TEST_CASE(unrestrictedFromOccOcc) {

    // Initialize a test set of orbital rotation generators for an occupied-occupied orbital space.
    GQCP::VectorX<double> kappa {3};
    kappa << 1, 2, 3;

    const GQCP::UOrbitalRotationGeneratorsComponent<double> occ_occ_generators_alpha {kappa};  // 3 (doubly-)occupied spatial orbitals.
    const GQCP::UOrbitalRotationGeneratorsComponent<double> occ_occ_generators_beta {kappa};   // 3 (doubly-)occupied spatial orbitals.

    // Construct the total generators in a larger orbital space.
    const auto full_generators_alpha = GQCP::UOrbitalRotationGeneratorsComponent<double>::FromOccOcc(occ_occ_generators_alpha, 4);  // 4 total spatial orbitals.
    const auto full_generators_beta = GQCP::UOrbitalRotationGeneratorsComponent<double>::FromOccOcc(occ_occ_generators_beta, 4);    // 4 total spatial orbitals.

    const GQCP::UOrbitalRotationGenerators<double> full_generators {full_generators_alpha, full_generators_beta};

    GQCP::SquareMatrix<double> full_kappa_matrix {4};
    // clang-format off
    full_kappa_matrix << 0, -1, -2, 0,
                         1,  0, -3, 0,
                         2,  3,  0, 0,
                         0,  0,  0, 0;
    // clang-format on
    BOOST_CHECK(full_generators.asMatrix().alpha().isApprox(full_kappa_matrix));
    BOOST_CHECK(full_generators.asMatrix().beta().isApprox(full_kappa_matrix));
}


/**
 *  Check if the `FromOccOcc` named constructor behaves correctly for the generalized case.
 */
BOOST_AUTO_TEST_CASE(generalizedFromOccOcc) {

    // Initialize a test set of orbital rotation generators for an occupied-occupied orbital space.
    GQCP::VectorX<double> kappa {3};
    kappa << 1, 2, 3;
    const GQCP::GOrbitalRotationGenerators<double> occ_occ_generators {kappa};  // 3 (doubly-)occupied spatial orbitals.

    // Construct the total generators in a larger orbital space.
    const auto full_generators = GQCP::GOrbitalRotationGenerators<double>::FromOccOcc(occ_occ_generators, 4);  // 4 total spatial orbitals.

    GQCP::SquareMatrix<double> full_kappa_matrix {4};
    // clang-format off
    full_kappa_matrix << 0, -1, -2, 0,
                         1,  0, -3, 0,
                         2,  3,  0, 0,
                         0,  0,  0, 0;
    // clang-format on
    BOOST_CHECK(full_generators.asMatrix().isApprox(full_kappa_matrix));
}
