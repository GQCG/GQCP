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

#define BOOST_TEST_MODULE "NuclearRepulsionOperator_test"

#include <boost/test/unit_test.hpp>

#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"


/**
 *  Check the internuclear repulsion energy for H2 from a reference value calculated by HORTON.
 */
BOOST_AUTO_TEST_CASE(NuclearRepulsion_h2) {

    const double ref_internuclear_repulsion_energy = 0.714285658963;

    // Create the dihydrogen nuclear framework.
    const auto molecule = GQCP::NuclearFramework::ReadXYZ("data/h2_szabo.xyz");
    BOOST_CHECK(std::abs(GQCP::NuclearRepulsionOperator(molecule).value() - ref_internuclear_repulsion_energy) < 1.0e-07);
}


/**
 *  Check the internuclear repulsion energy for H2O from a reference value calculated by HORTON.
 */
BOOST_AUTO_TEST_CASE(NuclearRepulsion_h2o) {

    const double ref_internuclear_repulsion_energy = 8.00236693455;

    // Create the water nuclear framework.
    const auto molecule = GQCP::NuclearFramework::ReadXYZ("data/h2o.xyz");
    BOOST_CHECK(std::abs(GQCP::NuclearRepulsionOperator(molecule).value() - ref_internuclear_repulsion_energy) < 1.0e-07);
}
