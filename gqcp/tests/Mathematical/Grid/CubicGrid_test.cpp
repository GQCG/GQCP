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

#define BOOST_TEST_MODULE "CubicGrid_test"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis/CartesianGTO.hpp"
#include "Mathematical/Grid/CubicGrid.hpp"
#include "Utilities/units.hpp"


/**
 *  Check if a cube file is made when writing field information.
 */
BOOST_AUTO_TEST_CASE(writeToCubeFile) {

    // Create a test molecule.
    const GQCP::Nucleus C {6, 0.0, 0.0, 0.0};
    const GQCP::Nucleus O {8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145)};  // from CCCBDB, STO-3G geometry
    const std::vector<GQCP::Nucleus> nuclei {C, O};
    const GQCP::Molecule molecule {{nuclei}};

    // Create a test scalar function to evaluate. In this case, we choose a GTO.
    const GQCP::CartesianExponents exponents {1, 0, 1};  // an x,z p-type GTO
    GQCP::Vector<double, 3> center;
    center << 5.0, 5.0, 5.5;
    const GQCP::CartesianGTO gto {1.0, exponents, center};

    // Set up a test cubic grid.
    GQCP::Vector<double, 3> origin = GQCP::Vector<double, 3>::Zero();
    const std::array<size_t, 3> steps {40, 40, 40};
    const std::array<double, 3> step_sizes {0.25, 0.25, 0.25};
    const GQCP::CubicGrid grid {origin, steps, step_sizes};

    // Evaluate the GTO on the cubic grid, and write the results to a cube file.
    const auto scalar_field = grid.evaluate(gto);
    grid.writeToCubeFile(scalar_field, "test_file.cube", molecule);
}


/**
 *  Check if reading in the data in an .rgrid-file is correct.
 */
BOOST_AUTO_TEST_CASE(ReadRegularGridFile) {

    // Read in the .rgrid-file, provide the reference values and check if it was parsed correctly.
    const auto grid = GQCP::CubicGrid::ReadRegularGridFile("data/benzene.rgrid");

    const GQCP::Vector<double, 3> ref_origin {-5.556445, -5.957031, -5.565039};
    const std::array<double, 3> ref_step_sizes {0.500586, 0.500586, 0.501367};
    const std::array<size_t, 3> ref_number_of_steps {24, 24, 24};


    // Check the results.
    BOOST_CHECK(grid.origin().isApprox(ref_origin, 1.0e-08));

    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(grid.numbersOfSteps(i) == ref_number_of_steps[i]);
        BOOST_CHECK(std::abs(grid.stepSize(i) - ref_step_sizes[i]) < 1.0e-08);
    }
}


/**
 *  Check if reading in the data in a GAUSSIAN Cube file is correct.
 */
BOOST_AUTO_TEST_CASE(ReadCubeFile) {

    // Read in the Cube, provide the reference values and check if it was parsed correctly.
    const auto grid = GQCP::CubicGrid::ReadCubeFile("data/benzene.cube");

    const GQCP::Vector<double, 3> ref_origin {-6.0, -6.0, -6.0};
    const std::array<double, 3> ref_step_sizes {0.2, 0.2, 0.2};
    const std::array<size_t, 3> ref_number_of_steps {60, 60, 60};


    // Check the results.
    BOOST_CHECK(grid.origin().isApprox(ref_origin, 1.0e-08));

    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(grid.numbersOfSteps(i) == ref_number_of_steps[i]);
        BOOST_CHECK(std::abs(grid.stepSize(i) - ref_step_sizes[i]) < 1.0e-08);
    }
}


/**
 *  Check if CubicGrid's integrate() method works as expected.
 */
BOOST_AUTO_TEST_CASE(integrate) {

    // Set up an example cubic grid.
    const GQCP::Vector<double, 3> origin = GQCP::Vector<double, 3>::Zero();
    const std::array<size_t, 3> number_of_steps {60, 60, 60};
    const std::array<double, 3> step_sizes {1.0, 1.0, 1.0};

    const GQCP::CubicGrid grid {origin, number_of_steps, step_sizes};


    // As an easy example, the scalar field f(r) = 1 should integrate to the total volume if the cubic grid.
    std::vector<double> field_values(grid.numberOfPoints(), 1.0);
    const GQCP::Field<double> scalar_field {field_values};

    BOOST_CHECK(std::abs(grid.integrate(scalar_field) - grid.totalVolume()) < 1.0e-12);
}


/**
 *  Check if an integration through CubicGrid equals an integration through an equivalent WeightedGrid.
 */



/**
 *  Check if the named constructor CubicGrid::Centered is correctly implemented
 */
BOOST_AUTO_TEST_CASE(Centered) {

    // We'll try to center a grid around (0,0,0) with 5 steps and a step size of 0.1.
    const size_t number_of_steps = 5;
    const double step_size = 0.1;
    const GQCP::Vector<double, 3> point = GQCP::Vector<double, 3>::Zero();

    const auto grid = GQCP::CubicGrid::Centered(point, number_of_steps, step_size);


    // Check if the origin of the grid appears to be correct.
    const GQCP::Vector<double, 3> ref_origin {-0.25, -0.25, -0.25};
    BOOST_CHECK(grid.origin().isApprox(ref_origin, 1.0e-08));
}
