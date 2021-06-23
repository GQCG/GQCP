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

#define BOOST_TEST_MODULE "elements"

#include <boost/test/unit_test.hpp>

#include "Molecule/elements.hpp"


/**
 *  Check a few conversions from element symbol to atomic number.
 */
BOOST_AUTO_TEST_CASE(elementToAtomicNumber) {

    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("H"), 1);
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("C"), 6);
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("Mg"), 12);
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("Nb"), 41);
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("Dy"), 66);
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("Pt"), 78);
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("U"), 92);
    BOOST_CHECK_EQUAL(GQCP::elements::elementToAtomicNumber("Og"), 118);
}


/**
 *  Check a few conversions from atomic number to element symbol.
 */
BOOST_AUTO_TEST_CASE(atomicNumberToElement) {

    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(1), "H");
    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(6), "C");
    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(12), "Mg");
    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(41), "Nb");
    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(66), "Dy");
    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(78), "Pt");
    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(92), "U");
    BOOST_CHECK_EQUAL(GQCP::elements::atomicNumberToElement(118), "Og");
}
