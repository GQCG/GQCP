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

#include "Basis/Transformations/ROrbitalRotationGenerators.hpp"
#include "Basis/Transformations/UOrbitalRotationGenerators.hpp"
#include "Basis/Transformations/UOrbitalRotationGeneratorsComponent.hpp"


/**
 *  Check the SimpleOrbitalRotationGenerators constructor from a kappa vector by using the ROrbitalRotationGenerators constructor.
 */
BOOST_AUTO_TEST_CASE(constructorVector) {

    // Check if the constructor throws upon receiving an argument of wrong dimensions.
    const GQCP::VectorX<double> kappa {4};  // '4' is not a triangular number.
    BOOST_CHECK_THROW(GQCP::ROrbitalRotationGenerators<double> generators {kappa}, std::invalid_argument);
}

/**
 *  Check the SimpleOrbitalRotationGenerators constructor from a kappa matrix by using the ROrbitalRotationGenerators constructor.
 */
BOOST_AUTO_TEST_CASE(constructorMatrix) {

    // Check if the constructor throws upon receiving an argument matrix that isn't anti-Hermitian.
    GQCP::SquareMatrix<double> kappa_matrix {3};
    // clang-format off
    kappa_matrix << 0,  1,  2,
                    1,  0,  3,
                    2,  3,  0;
    // clang-format on
    BOOST_CHECK_THROW(GQCP::ROrbitalRotationGenerators<double> generators {kappa_matrix}, std::invalid_argument);
}


/**
 *  Check if the `asMatrix` API correctly converts ROrbitalRotationGenerators to their matrix representation.
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
    const GQCP::ROrbitalRotationGenerators<double> generators {kappa};
    BOOST_CHECK(generators.asMatrix().isApprox(kappa_matrix));
}


/**
 *  Check if the `FromOccupationType` named constructor behaves correctly for the restricted case.
 */
BOOST_AUTO_TEST_CASE(FromOccupationType) {

    // Initialize a test set of orbital rotation generators for an occupied-occupied orbital space.
    GQCP::VectorX<double> kappa {3};
    kappa << 1, 2, 3;
    const GQCP::ROrbitalRotationGenerators<double> occ_occ_generators {kappa};  // 3 (doubly-)occupied spatial orbitals.

    // Construct the total generators in a larger orbital space.
    const auto full_generators = GQCP::ROrbitalRotationGenerators<double>::FromOccupationTypes(occ_occ_generators, GQCP::OccupationType::k_occupied, GQCP::OccupationType::k_occupied, 4);  // 4 total spatial orbitals.

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
 *  Check if the `FromOccupationType` named constructor behaves correctly for the restricted case.
 */
// BOOST_AUTO_TEST_CASE(FromOccupationTypeOccupiedVirtual) {

//     // Initialize a test set of orbital rotation generators for an occupied-occupied orbital space.
//     GQCP::VectorX<double> kappa {3};
//     kappa << 1, 2, 3;
//     const GQCP::ROrbitalRotationGenerators<double> occ_vir_generators {kappa};  // 3 (doubly-)occupied spatial orbitals.

//     // Construct the total generators in a larger orbital space.
//     const auto full_generators = GQCP::ROrbitalRotationGenerators<double>::FromOccupationTypes(occ_vir_generators, GQCP::OccupationType::k_occupied, GQCP::OccupationType::k_virtual, 5);  // 5 total spatial orbitals.

//     GQCP::SquareMatrix<double> full_kappa_matrix {5};
//     // clang-format off
//     full_kappa_matrix << 0, 0, 0, -1, -2,
//                          0, 0, 1,  0, -3,
//                          0, 0, 2,  3,  0,
//                          0, 0, 0,  0,  0,
//                          0, 0, 0,  0,  0;
//     // clang-format on
//     BOOST_CHECK(full_generators.asMatrix().isApprox(full_kappa_matrix));
// }


/**
 *  Check if the `FromOccupationType` named constructor behaves correctly for the unrestricted case.
 */
BOOST_AUTO_TEST_CASE(FromOccupationTypeUnrestricted) {

    // Initialize a test set of orbital rotation generators for an occupied-occupied orbital space.
    GQCP::VectorX<double> kappa {3};
    kappa << 1, 2, 3;

    const GQCP::UOrbitalRotationGeneratorsComponent<double> occ_occ_generators_alpha {kappa};  // 3 (doubly-)occupied spin-orbitals.
    const GQCP::UOrbitalRotationGeneratorsComponent<double> occ_occ_generators_beta {kappa};   // 3 (doubly-)occupied spin-orbitals.

    const GQCP::UOrbitalRotationGenerators<double> occ_occ_generators {occ_occ_generators_alpha, occ_occ_generators_beta};

    // Construct the total generators in a larger orbital space.
    const auto full_generators = GQCP::UOrbitalRotationGenerators<double>::FromOccupationTypes(occ_occ_generators, GQCP::OccupationType::k_occupied, GQCP::OccupationType::k_occupied, 4);  // 4 total spin-orbitals.

    GQCP::SquareMatrix<double> full_kappa_matrix {4};
    // clang-format off
    full_kappa_matrix << 0, -1, -2, 0,
                         1,  0, -3, 0,
                         2,  3,  0, 0,
                         0,  0,  0, 0;
    // clang-format on
    BOOST_CHECK(full_generators.alpha().asMatrix().isApprox(full_kappa_matrix));
    BOOST_CHECK(full_generators.beta().asMatrix().isApprox(full_kappa_matrix));
}
