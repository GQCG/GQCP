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

#define BOOST_TEST_MODULE "Nucleus"

#include <boost/test/unit_test.hpp>

#include "Molecule/Nucleus.hpp"


BOOST_AUTO_TEST_CASE(Nucleus_constructor) {

    GQCP::Nucleus nucleus {1, 0.0, 0.1, 0.2};
}


BOOST_AUTO_TEST_CASE(Nucleus_sortComparer) {

    GQCP::Nucleus nucleus1 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus2 {2, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus3 {2, 0.1, 0.2, 0.2};


    // Check if sortComparer does what is expected
    const auto sort_comparer = GQCP::Nucleus::sortComparer();
    BOOST_CHECK(sort_comparer(nucleus1, nucleus3));  // x1 < x2
    BOOST_CHECK(sort_comparer(nucleus2, nucleus3));  // x1 < x2

    BOOST_CHECK(!(sort_comparer(nucleus1, nucleus2)));  // same positions, result should be false

    // Check if the tolerance works
    const auto sort_comparer_tolerance = GQCP::Nucleus::sortComparer(0.2);
    BOOST_CHECK(!(sort_comparer_tolerance(nucleus2, nucleus3)));
    BOOST_CHECK(!(sort_comparer_tolerance(nucleus3, nucleus2)));
}


BOOST_AUTO_TEST_CASE(Nucleus_equalityComparer) {

    GQCP::Nucleus nucleus1 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus2 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus3 {2, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus4 {1, 0.1, 0.2, 0.3};

    // Check if equalityComparer does what is expected
    const auto equality_comparer = GQCP::Nucleus::equalityComparer();
    BOOST_CHECK(equality_comparer(nucleus1, nucleus2));

    // Check if different charges cause inequality
    BOOST_CHECK(!(equality_comparer(nucleus1, nucleus3)));

    // Check if different coordinates cause inequality
    BOOST_CHECK(!(equality_comparer(nucleus1, nucleus4)));


    // Check if the tolerance works
    const auto equality_comparer_tolerance = GQCP::Nucleus::equalityComparer(0.2);
    BOOST_CHECK(equality_comparer_tolerance(nucleus1, nucleus2));
}


BOOST_AUTO_TEST_CASE(Nucleus_operator_ostream) {

    GQCP::Nucleus nucleus1 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus2 {2, 0.1, 0.2, 0.3};


    std::cout << nucleus1 << std::endl;
    std::cout << nucleus2 << std::endl;
}


BOOST_AUTO_TEST_CASE(calculateDistanceWith) {

    // Create some nuclei
    GQCP::Nucleus nucleus1 {1, 0, 3, 0};
    GQCP::Nucleus nucleus2 {1, 0, 0, 4};
    GQCP::Nucleus nucleus3 {1, 3, 0, 0};
    GQCP::Nucleus nucleus4 {1, 0, 0, 5};


    // Check their distances
    BOOST_CHECK(std::abs(nucleus1.calculateDistanceWith(nucleus2) - 5) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus1.calculateDistanceWith(nucleus3) - std::sqrt(18.0)) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus1.calculateDistanceWith(nucleus2) - nucleus2.calculateDistanceWith(nucleus3)) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus2.calculateDistanceWith(nucleus3) - 5) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus2.calculateDistanceWith(nucleus4) - 1) < 1.0e-12);

    // Check that the distances are symmetric
    BOOST_CHECK(std::abs(nucleus1.calculateDistanceWith(nucleus2) - nucleus2.calculateDistanceWith(nucleus1)) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus1.calculateDistanceWith(nucleus3) - nucleus3.calculateDistanceWith(nucleus1)) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus2.calculateDistanceWith(nucleus3) - nucleus3.calculateDistanceWith(nucleus2)) < 1.0e-12);
}
