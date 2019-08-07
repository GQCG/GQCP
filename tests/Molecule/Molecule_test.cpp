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
#define BOOST_TEST_MODULE "Molecule"

#include <boost/test/included/unit_test.hpp>

#include "Molecule/Molecule.hpp"


BOOST_AUTO_TEST_CASE ( constructor_nuclei_charge ) {

    // Create a fictitious molecule from some nuclei (charge, x, y ,z)
    std::vector<GQCP::Nucleus> nuclei = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    // Check if we can create any anion
    GQCP::Molecule molecule2 (nuclei, -2);


    // Check if we can't create a cation with charge larger than the nucleic charge
    BOOST_CHECK_NO_THROW(GQCP::Molecule (nuclei, +3));
    BOOST_CHECK_THROW(GQCP::Molecule (nuclei, +11), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( molecule_ion_constructor ) {

    // Create some Molecule objects
    const std::string xyzfilename = "data/h2o.xyz";
    auto water = GQCP::Molecule::ReadXYZ(xyzfilename);
    auto water_anion = GQCP::Molecule::ReadXYZ(xyzfilename, -1);
    auto water_neutral = GQCP::Molecule::ReadXYZ(xyzfilename, 0);
    auto water_cation = GQCP::Molecule::ReadXYZ(xyzfilename, +1);

    // Test the number of electrons created by the constructor
    BOOST_CHECK_EQUAL(water.numberOfElectrons(), 10);
    BOOST_CHECK_EQUAL(water_anion.numberOfElectrons(), 11);
    BOOST_CHECK_EQUAL(water_neutral.numberOfElectrons(), 10);
    BOOST_CHECK_EQUAL(water_cation.numberOfElectrons(), 9);
}


BOOST_AUTO_TEST_CASE ( Molecule_operator_ostream ) {

    std::vector<GQCP::Nucleus> nuclei = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };
    GQCP::Molecule molecule (nuclei);

    std::cout << molecule << std::endl;
}


BOOST_AUTO_TEST_CASE ( HChain ) {

    GQCP::Molecule h_chain = GQCP::Molecule::HChain(3, 1.0);
    BOOST_CHECK(h_chain.numberOfAtoms() == 3);
    BOOST_CHECK(h_chain.numberOfElectrons() == 3);

    GQCP::Molecule h_chain_charged = GQCP::Molecule::HChain(4, 1.5, +2);
    BOOST_CHECK(h_chain_charged.numberOfAtoms() == 4);
    BOOST_CHECK(h_chain_charged.numberOfElectrons() == 2);
}


BOOST_AUTO_TEST_CASE ( H2Chain ) {

    GQCP::Molecule h2_chain = GQCP::Molecule::H2Chain(2, 1.0, 1.5);
    BOOST_CHECK(h2_chain.numberOfAtoms() == 4);
    BOOST_CHECK(h2_chain.numberOfElectrons() == 4);

    GQCP::Molecule h2_chain_charged = GQCP::Molecule::H2Chain(4, 2.0, 3.8, -2);
    BOOST_CHECK(h2_chain_charged.numberOfAtoms() == 8);
    BOOST_CHECK(h2_chain_charged.numberOfElectrons() == 10);
}
