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

#define BOOST_TEST_MODULE "units"

#include <boost/test/unit_test.hpp>

#include "Utilities/units.hpp"


/**
 *  Check if the conversion from Bohr to Angstrom is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(bohr_angstrom) {

    // 1 bohr = 0.529 Angstrom
    BOOST_CHECK(std::abs((GQCP::units::bohr_to_angstrom(1) - 0.529)) < 1.0e-03);

    // Check if back-to-back transformations have no effect.
    BOOST_CHECK(std::abs((GQCP::units::bohr_to_angstrom(GQCP::units::angstrom_to_bohr(0.5)) - 0.5)) < 1.0e-12);
}
