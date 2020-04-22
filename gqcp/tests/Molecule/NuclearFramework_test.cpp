// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2019  the GQCG developers
//
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
//
#define BOOST_TEST_MODULE "NuclearFramework"

#include <boost/test/unit_test.hpp>

#include "Molecule/NuclearFramework.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Test if we can't create a NuclearFramework with duplicate nuclei.
 */
BOOST_AUTO_TEST_CASE(duplicate_nuclei_constructor) {

    // Make some nuclei
    const GQCP::Nucleus nucleus1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus nucleus2 {1, 1.0, 0.0, 0.0};

    const std::vector<GQCP::Nucleus> nuclei1 {nucleus1, nucleus1};
    const std::vector<GQCP::Nucleus> nuclei2 {nucleus1, nucleus2};


    // Check if we can't create a NuclearFramework with duplicate nuclei
    BOOST_CHECK_THROW(GQCP::NuclearFramework nuclear_framework {nuclei1}, std::invalid_argument);

    // Check if a correct argument doesn't throw
    BOOST_CHECK_NO_THROW(GQCP::NuclearFramework nuclear_framework {nuclei2});
}


/**
 *  Check if totalNucleicCharge() works as expected.
 */
BOOST_AUTO_TEST_CASE(totalNucleicCharge) {

    // Create a fictitious molecule from some nuclei (charge, x, y ,z)
    const std::vector<GQCP::Nucleus> nuclei = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}};

    const GQCP::NuclearFramework nuclear_framework {nuclei};
    BOOST_CHECK_EQUAL(nuclear_framework.totalNucleicCharge(), 10);
}


/**
 *  Check if internuclearDistance() works as expected.
 */
BOOST_AUTO_TEST_CASE(internuclearDistance) {

    // Create a fictitious molecule from some nuclei (charge, x, y ,z)
    const std::vector<GQCP::Nucleus> nuclei = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}};
    const GQCP::NuclearFramework nuclear_framework {nuclei};


    // Check if we get throws when the indices are out of bounds
    BOOST_CHECK_THROW(nuclear_framework.internuclearDistance(0, 5), std::invalid_argument);
    BOOST_CHECK_THROW(nuclear_framework.internuclearDistance(8, 2), std::invalid_argument);

    // Check if we don't get throws when the indices behave correctly
    BOOST_CHECK_NO_THROW(nuclear_framework.internuclearDistance(0, 0));

    // Check if the function works
    BOOST_CHECK(std::abs(nuclear_framework.internuclearDistance(1, 3) - 1) < 1.0e-12);
}


/**
 *  Check if the basic methods work on H2.
 */
BOOST_AUTO_TEST_CASE(methods_h2) {

    // Create the dihydrogen nuclear framework
    const auto h2 = GQCP::NuclearFramework::ReadXYZ("data/h2_szabo.xyz");

    // Test the basic methods
    BOOST_CHECK_EQUAL(h2.numberOfNuclei(), 2);
    BOOST_CHECK_EQUAL(h2.totalNucleicCharge(), 2);
}


/**
 *  Check if the basic methods work on H2O.
 */
BOOST_AUTO_TEST_CASE(methods_water) {

    // Create the water nuclear framework
    const auto water = GQCP::NuclearFramework::ReadXYZ("data/h2o.xyz");

    // Test the basic methods
    BOOST_CHECK_EQUAL(water.numberOfNuclei(), 3);
    BOOST_CHECK_EQUAL(water.totalNucleicCharge(), 10);
}


/**
 *  Check if the NuclearFramework::HChain throws when expected.
 */
BOOST_AUTO_TEST_CASE(HChain_throws) {

    BOOST_CHECK_THROW(GQCP::NuclearFramework::HChain(0, 1.0), std::invalid_argument);   // can't create 0 H-nuclei
    BOOST_CHECK_THROW(GQCP::NuclearFramework::HChain(1, -1.0), std::invalid_argument);  // can't have negative spacing
}


/**
 *  Check if the NuclearFramework::H2Chain throws when expected.
 */
BOOST_AUTO_TEST_CASE(H2Chain_throws) {

    BOOST_CHECK_THROW(GQCP::NuclearFramework::H2Chain(0, 1.0, 2.0), std::invalid_argument);   // can't create 0 H2-molecules
    BOOST_CHECK_THROW(GQCP::NuclearFramework::H2Chain(1, -1.0, 1.0), std::invalid_argument);  // can't have negative spacing
    BOOST_CHECK_THROW(GQCP::NuclearFramework::H2Chain(1, 1.0, -1.0), std::invalid_argument);  // can't have negative spacing
}


/**
 *  Check if the construction of an H-chain works, with a manual example.
 */
BOOST_AUTO_TEST_CASE(HChain) {

    // Check the construction for a H3-chain.
    const GQCP::NuclearFramework h_chain1 = GQCP::NuclearFramework::HChain(3, 1.0);
    BOOST_CHECK(h_chain1.numberOfNuclei() == 3);
    BOOST_CHECK(h_chain1.totalNucleicCharge() == 3);

    BOOST_CHECK(std::abs(h_chain1.internuclearDistance(0, 1) - 1.0) < 1.0e-12);
    BOOST_CHECK(std::abs(h_chain1.internuclearDistance(0, 2) - 2.0) < 1.0e-12);


    // Check the construction for a H4-chain.
    const GQCP::NuclearFramework h_chain2 = GQCP::NuclearFramework::HChain(4, 1.5);
    BOOST_CHECK(h_chain1.numberOfNuclei() == 3);
    BOOST_CHECK(h_chain2.totalNucleicCharge() == 4);

    BOOST_CHECK(std::abs(h_chain2.internuclearDistance(0, 1) - 1.5) < 1.0e-12);
    BOOST_CHECK(std::abs(h_chain2.internuclearDistance(0, 2) - 3.0) < 1.0e-12);
    BOOST_CHECK(std::abs(h_chain2.internuclearDistance(0, 3) - 4.5) < 1.0e-12);
    BOOST_CHECK(std::abs(h_chain2.internuclearDistance(1, 2) - 1.5) < 1.0e-12);
    BOOST_CHECK(std::abs(h_chain2.internuclearDistance(1, 3) - 3.0) < 1.0e-12);
    BOOST_CHECK(std::abs(h_chain2.internuclearDistance(2, 3) - 1.5) < 1.0e-12);
}


/**
 *  Check if the construction of an H2-chain works, with a manual example.
 */
BOOST_AUTO_TEST_CASE(H2Chain) {

    // Check the construction for 2 H2-molecules
    const GQCP::NuclearFramework h2_chain1 = GQCP::NuclearFramework::H2Chain(2, 1.0, 1.5);
    BOOST_CHECK(h2_chain1.numberOfNuclei() == 4);
    BOOST_CHECK(h2_chain1.totalNucleicCharge() == 4);

    BOOST_CHECK(std::abs(h2_chain1.internuclearDistance(0, 1) - 1.0) < 1.0e-12);
    BOOST_CHECK(std::abs(h2_chain1.internuclearDistance(0, 2) - 2.5) < 1.0e-12);
    BOOST_CHECK(std::abs(h2_chain1.internuclearDistance(0, 3) - 3.5) < 1.0e-12);
    BOOST_CHECK(std::abs(h2_chain1.internuclearDistance(1, 2) - 1.5) < 1.0e-12);
    BOOST_CHECK(std::abs(h2_chain1.internuclearDistance(1, 3) - 2.5) < 1.0e-12);
    BOOST_CHECK(std::abs(h2_chain1.internuclearDistance(2, 3) - 1.0) < 1.0e-12);


    // Check the construction for 4 H2-molecules
    const GQCP::NuclearFramework h2_chain2 = GQCP::NuclearFramework::H2Chain(4, 2.0, 3.8);
    BOOST_CHECK(h2_chain2.numberOfNuclei() == 8);
    BOOST_CHECK(h2_chain2.totalNucleicCharge() == 8);
}


/**
 *  Check if the NuclearFramework::HRingFromRadius throws when expected.
 */
BOOST_AUTO_TEST_CASE(HRingFromRadius_throws) {

    BOOST_CHECK_THROW(GQCP::NuclearFramework::HRingFromRadius(0, 1.0), std::invalid_argument);   // can't create 0 H-nuclei
    BOOST_CHECK_THROW(GQCP::NuclearFramework::HRingFromRadius(1, -1.0), std::invalid_argument);  // can't a negative radius
}


/**
 *  Check if the construction of an H4-ring from a given circumscribing radius works as expected.
 */
BOOST_AUTO_TEST_CASE(HRingFromRadius) {

    const auto ring = GQCP::NuclearFramework::HRingFromRadius(4, 1.0);
    const auto distance = boost::math::constants::root_two<double>();  // 2 R sin(pi/n)

    BOOST_CHECK(ring.numberOfNuclei() == 4);
    BOOST_CHECK(ring.totalNucleicCharge() == 4);

    BOOST_CHECK(std::abs(ring.internuclearDistance(0, 1) - distance) < 1.0e-12);
    BOOST_CHECK(std::abs(ring.internuclearDistance(1, 2) - distance) < 1.0e-12);
    BOOST_CHECK(std::abs(ring.internuclearDistance(2, 3) - distance) < 1.0e-12);
    BOOST_CHECK(std::abs(ring.internuclearDistance(3, 0) - distance) < 1.0e-12);
}


/**
 *  Check if the construction of an H4-ring from a given neighbour distance works as expected.
 */
BOOST_AUTO_TEST_CASE(HRingFromDistance) {

    const double distance = 1.0;
    const auto ring = GQCP::NuclearFramework::HRingFromDistance(4, distance);

    BOOST_CHECK(ring.numberOfNuclei() == 4);
    BOOST_CHECK(ring.totalNucleicCharge() == 4);

    BOOST_CHECK(std::abs(ring.internuclearDistance(0, 1) - distance) < 1.0e-12);
    BOOST_CHECK(std::abs(ring.internuclearDistance(1, 2) - distance) < 1.0e-12);
    BOOST_CHECK(std::abs(ring.internuclearDistance(2, 3) - distance) < 1.0e-12);
    BOOST_CHECK(std::abs(ring.internuclearDistance(3, 0) - distance) < 1.0e-12);
}
