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

#define BOOST_TEST_MODULE "LondonCartesianGTO"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Functions/LondonCartesianGTO.hpp"


/**
 *  Check if the London phase factor is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(phase_factor) {

    // Set up toy vectors.
    const GQCP::Vector<double, 3> K {1.0, 0.0, -1.0};  // The center.
    const GQCP::Vector<double, 3> G {1.0, 2.0, 3.0};   // The gauge origin.
    const GQCP::Vector<double, 3> B {1.0, 1.0, 1.0};   // The magnetic field.

    const GQCP::Vector<double, 3> r {-1.0, 1.0, -1.0};  // The position at which to evaluate.


    // Initialize a `LondonCartesianGTO` and check if its `phaseFactor` is correctly implemented.
    const GQCP::HomogeneousMagneticField field {B, G};

    const GQCP::CartesianExponents exponents {0, 0, 0};
    const GQCP::CartesianGTO gto {1.0, exponents, K};
    const GQCP::LondonCartesianGTO london_gto {field, gto};

    std::complex<double> ref {-0.653643620863611, 0.756802495307928};  // This is equal to exp(-4i).
    BOOST_CHECK(std::abs(london_gto.phaseFactor(r) - ref) < 1.0e-14);
}
