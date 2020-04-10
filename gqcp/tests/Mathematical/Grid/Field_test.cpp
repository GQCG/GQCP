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
#define BOOST_TEST_MODULE "EvaluatedScalarFuctionOnCubicGrid"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


#include "Grid/EvaluatedScalarFunctionOnCubicGrid.hpp"
#include "units.hpp"
#include "Basis/CartesianGTO.hpp"


BOOST_AUTO_TEST_CASE ( EvaluatedScalarFunctionOnCubicGrid_test_placeholder ) {

    GQCP::Atom C (6, 0.0, 0.0, 0.0);
    GQCP::Atom O (8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Atom> atoms {C, O};
    GQCP::Molecule CO (atoms);

    GQCP::CartesianExponents exponents1 (1, 0, 1);
    GQCP::Vector<double, 3> center1;
    center1 << 5.0, 5.0, 5.5;
    GQCP::Vector<double, 3> r1;
    r1 << 0.0, 1.0, 0.0;
    GQCP::CartesianGTO gto1 (1.0, exponents1, center1);

    GQCP::CubicGrid grid (GQCP::Vector<double,3>::Zero(3), {40, 40, 40}, {0.25, 0.25, 0.25});

    GQCP::EvaluatedScalarFunctionOnCubicGrid gto1_evaluated_on_cubic_grid(grid, gto1);
    gto1_evaluated_on_cubic_grid.toCubeFile("gqcp_fun_times.cube", CO);

    BOOST_CHECK(true);
}
