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
#define BOOST_TEST_MODULE "NuclearRepulsionOperator_test"

#include <boost/test/unit_test.hpp>

#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"


BOOST_AUTO_TEST_CASE ( NuclearRepulsion_h2 ) {

    // We have a reference internuclear repulsion energy from HORTON
    double ref_internuclear_repulsion_energy = 0.714285658963;

    // Create the dihydrogen nuclear framework
    auto h2 = GQCP::NuclearFramework::ReadXYZ("data/h2_szabo.xyz");

    // Test the calculation of the nuclear repulsion energy for hydrogen 
    BOOST_CHECK(std::abs(GQCP::NuclearRepulsionOperator(h2).value() - ref_internuclear_repulsion_energy) < 1.0e-07);  // reference data from horton
}


BOOST_AUTO_TEST_CASE ( NuclearRepulsion_h2o ) {

    // We have reference internuclear repulsion energy from HORTON
    double ref_internuclear_repulsion_energy = 8.00236693455;

    // Create the water nuclear framework
    auto water = GQCP::NuclearFramework::ReadXYZ("data/h2o.xyz");

    // Test the calculation of the nuclear repulsion energy
    BOOST_CHECK(std::abs(GQCP::NuclearRepulsionOperator(water).value() - ref_internuclear_repulsion_energy) < 1.0e-07);  // reference data from horton
}
