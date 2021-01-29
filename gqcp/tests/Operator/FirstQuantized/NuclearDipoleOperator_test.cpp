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

#define BOOST_TEST_MODULE "NuclearDipoleOperator_test"

#include <boost/test/unit_test.hpp>

#include "Operator/FirstQuantized/NuclearDipoleOperator.hpp"


/**
 *  Check the nuclear dipole moment for a toy nuclear framework.
 */
BOOST_AUTO_TEST_CASE(NuclearDipole) {

    const GQCP::Nucleus H {1, 0, 1, 2};
    const GQCP::Nucleus O {8, 2, 4, 8};
    const GQCP::NuclearFramework nuclear_framework {{H, O}};

    const auto dipole_moment = GQCP::NuclearDipoleOperator(nuclear_framework).value();
    BOOST_CHECK(dipole_moment.isApprox(GQCP::Vector<double, 3> {16, 33, 66}));
}
