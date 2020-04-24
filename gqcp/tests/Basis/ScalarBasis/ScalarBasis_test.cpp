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

#define BOOST_TEST_MODULE "ScalarBasis"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/Operator.hpp"


BOOST_AUTO_TEST_CASE(Scalar_basis_constructor) {

    // Check if we can construct an ScalarBasis<GTOShell> object
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::ScalarBasis<GQCP::GTOShell> basis {water, "STO-3G"};
}


BOOST_AUTO_TEST_CASE(numberOfBasisFunctions) {

    // Check the number of basis functions in water
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::ScalarBasis<GQCP::GTOShell> basis {water, "STO-3G"};

    BOOST_CHECK_EQUAL(basis.numberOfBasisFunctions(), 7);
}


/**
 *  Check if the underlying shell set created from a molecule is the same as the one created from a nuclear framework
 */
BOOST_AUTO_TEST_CASE(molecule_nuclear_framework) {

    // Create the two scalar bases and check if the underlying shell sets are equal
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");

    GQCP::ScalarBasis<GQCP::GTOShell> basis_molecule {water, "STO-3G"};
    GQCP::ScalarBasis<GQCP::GTOShell> basis_nuclear_framework {water.nuclearFramework(), "STO-3G"};

    BOOST_CHECK(basis_molecule.shellSet().asVector() == basis_nuclear_framework.shellSet().asVector());
}


BOOST_AUTO_TEST_CASE(dissociatedMoleculeBasis) {

    // Test if we can succesfully initialize NO+ at long intra molecular distance
    auto N = GQCP::Nucleus(7, 3.5, 0, 0);
    auto O = GQCP::Nucleus(8, -3.5, 0, 0);
    std::vector<GQCP::Nucleus> nuclei {N, O};
    auto NO = GQCP::Molecule(nuclei, +1);
    GQCP::ScalarBasis<GQCP::GTOShell> basis {NO, "STO-3G"};

    BOOST_CHECK_NO_THROW(GQCP::ScalarBasis<GQCP::GTOShell> basis (NO, "STO-3G"));
}
