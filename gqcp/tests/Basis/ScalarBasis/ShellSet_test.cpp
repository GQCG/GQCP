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

#define BOOST_TEST_MODULE "ShellSet"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"


/**
 *  Check if basic methods work as expected.
 */
BOOST_AUTO_TEST_CASE(basic) {

    // Generate a scalar basis of GTOs with at least a p-type function inside.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const auto shellset = GQCP::GTOBasisSet("STO-3G").generate(molecule);


    /**
     *  H2O in an STO-3G basisset contains:
     *      - 5 shells (sets of functions on the same nucleuswith the same angular momentum): 1s on H (2 times), 1s on O, 2s on O, 3s on O
     *      - 7 basis functions: 1 1s on H (2 times), 1 1s on O, 1 2s on O and 3 2p on O
     */
    BOOST_CHECK_EQUAL(shellset.numberOfShells(), 5);
    BOOST_CHECK_EQUAL(shellset.numberOfBasisFunctions(), 7);
    BOOST_CHECK_EQUAL(shellset.maximumNumberOfPrimitives(), 3);  // 3 primitives for O's p-type GTO
    BOOST_CHECK_EQUAL(shellset.maximumAngularMomentum(), 1);     // O has a p-type basis function
}
