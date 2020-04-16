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
#define BOOST_TEST_MODULE "Field"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain

#include "Basis/ScalarBasis/CartesianGTO.hpp"
#include "Mathematical/Grid/CubicGrid.hpp"
#include "Mathematical/Grid/Field.hpp"
#include "Utilities/units.hpp"


/**
 *  Check if a cube file is made when writing field information.
 */
BOOST_AUTO_TEST_CASE ( toCubeFile ) {

    // Create a test molecule.
    const GQCP::Nucleus C (6, 0.0, 0.0, 0.0);
    const GQCP::Nucleus O (8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145));  // from CCCBDB, STO-3G geometry
    const std::vector<GQCP::Nucleus> nuclei {C, O};
    const GQCP::Molecule molecule (nuclei);

    // Create a test scalar function to evaluate. In this case, we choose a GTO.
    const GQCP::CartesianExponents exponents (1, 0, 1);  // an x,z p-type GTO
    GQCP::Vector<double, 3> center;
    center << 5.0, 5.0, 5.5;
    const GQCP::CartesianGTO gto (1.0, exponents, center);

    // Set up a test cubic grid.
    GQCP::Vector<double, 3> origin = GQCP::Vector<double, 3>::Zero();
    const std::array<size_t, 3> steps = {40, 40, 40};
    const std::array<double, 3> step_sizes = {0.25, 0.25, 0.25};
    GQCP::CubicGrid grid (origin, steps, step_sizes);

    // Evaluate the GTO on the cubic grid, and write the results to a cube file.
    grid.evaluate(gto).toCubeFile("test_file.cube", molecule);
}
