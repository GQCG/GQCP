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

#define BOOST_TEST_MODULE "Molecule"

#include <boost/test/unit_test.hpp>

#include "Molecule/Molecule.hpp"


/**
 *  Check if constructors related to charged species behave as expected.
 */
BOOST_AUTO_TEST_CASE(constructor_nuclei_charge) {

    // Create a fictitious molecule from some nuclei (charge, x, y ,z)
    const std::vector<GQCP::Nucleus> nuclei = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}};

    // Check if we can create any anion
    const GQCP::Molecule molecule2(nuclei, -2);


    // Check if we can't create a cation with charge larger than the nucleic charge
    BOOST_CHECK_NO_THROW(GQCP::Molecule(nuclei, +3));
    BOOST_CHECK_THROW(GQCP::Molecule(nuclei, +11), std::invalid_argument);
}


/**
 *  Check if numberOfElectrons() is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(molecule_ion_constructor) {

    // Create some Molecule objects
    const std::string xyzfilename = "data/h2o.xyz";
    const auto water = GQCP::Molecule::ReadXYZ(xyzfilename);
    const auto water_anion = GQCP::Molecule::ReadXYZ(xyzfilename, -1);
    const auto water_neutral = GQCP::Molecule::ReadXYZ(xyzfilename, 0);
    const auto water_cation = GQCP::Molecule::ReadXYZ(xyzfilename, +1);

    // Test the number of electrons created by the constructor
    BOOST_CHECK_EQUAL(water.numberOfElectrons(), 10);
    BOOST_CHECK_EQUAL(water_anion.numberOfElectrons(), 11);
    BOOST_CHECK_EQUAL(water_neutral.numberOfElectrons(), 10);
    BOOST_CHECK_EQUAL(water_cation.numberOfElectrons(), 9);
}


/**
 *  Check if the API still supports Molecule<<
 */
BOOST_AUTO_TEST_CASE(Molecule_operator_ostream) {

    const std::vector<GQCP::Nucleus> nuclei = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}};
    const GQCP::Molecule molecule {nuclei};

    std::cout << molecule << std::endl;
}


/**
 *  Check if the correct number of electrons is implemented for H-chains.
 */
BOOST_AUTO_TEST_CASE(HChain) {

    const auto h_chain = GQCP::Molecule::HChain(3, 1.0);
    BOOST_CHECK(h_chain.numberOfAtoms() == 3);
    BOOST_CHECK(h_chain.numberOfElectrons() == 3);

    const auto h_chain_charged = GQCP::Molecule::HChain(4, 1.5, +2);
    BOOST_CHECK(h_chain_charged.numberOfAtoms() == 4);
    BOOST_CHECK(h_chain_charged.numberOfElectrons() == 2);
}


/**
 *  Check if the correct number of electrons is implemented for H2-chains.
 */
BOOST_AUTO_TEST_CASE(H2Chain) {

    const auto h2_chain = GQCP::Molecule::H2Chain(2, 1.0, 1.5);
    BOOST_CHECK(h2_chain.numberOfAtoms() == 4);
    BOOST_CHECK(h2_chain.numberOfElectrons() == 4);

    const auto h2_chain_charged = GQCP::Molecule::H2Chain(4, 2.0, 3.8, -2);
    BOOST_CHECK(h2_chain_charged.numberOfAtoms() == 8);
    BOOST_CHECK(h2_chain_charged.numberOfElectrons() == 10);
}


/**
 *  Check if the molecular charge and related methods are correctly implemented.
 */
BOOST_AUTO_TEST_CASE(charge) {

    const auto h_chain = GQCP::Molecule::HChain(3, 1.0);  // neutral species
    BOOST_CHECK(h_chain.totalNucleicCharge() == 3);
    BOOST_CHECK(h_chain.numberOfElectrons() == 3);
    BOOST_CHECK(h_chain.numberOfElectronPairs() == 1);  // should round down
    BOOST_CHECK(h_chain.charge() == 0);


    const auto h_chain_cation = GQCP::Molecule::HChain(3, 1.0, +2);  // charge +2
    BOOST_CHECK(h_chain_cation.totalNucleicCharge() == 3);
    BOOST_CHECK(h_chain_cation.numberOfElectrons() == 1);
    BOOST_CHECK(h_chain_cation.numberOfElectronPairs() == 0);  // should round down
    BOOST_CHECK(h_chain_cation.charge() == +2);

    const auto h_chain_anion = GQCP::Molecule::HChain(3, 1.0, -2);  // charge -2
    BOOST_CHECK(h_chain_anion.totalNucleicCharge() == 3);
    BOOST_CHECK(h_chain_anion.numberOfElectrons() == 5);
    BOOST_CHECK(h_chain_anion.numberOfElectronPairs() == 2);  // should round down
    BOOST_CHECK(h_chain_anion.charge() == -2);
}
