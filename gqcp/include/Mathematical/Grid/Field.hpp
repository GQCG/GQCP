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

#pragma once


#include "Molecule/Molecule.hpp"
#include "Utilities/type_traits.hpp"

#include <string>
#include <vector>


namespace GQCP {


/**
 *  A set of function values corresponding to points in space.
 * 
 *  @tparam T           the type of the evaluated function values
 *  @tparam _Grid       the type of grid that is used to represent the points in space
 */
template <typename T, typename _Grid>
class Field {
public:
    using Grid = _Grid;


private:
    std::vector<T> values;  // the evaluated function values, in the order of the grid's loop
    Grid grid;              // the grid associated with this field


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param values           the evaluated function values, in the order of the grid's loop
     *  @param grid             the grid associated with this field
     */
    Field(const std::vector<T>& values, const Grid& grid) :
        values {values},
        grid {grid} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Write this field's values to a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/).
     *
     *  @param filename     the name of the cubefile that has to be generated
     *  @param molecule     the molecule that should be placed in the cubefile
     */
    void toCubeFile(const std::string& filename, const Molecule& molecule) const {

        // Prepare some variables.
        std::ofstream cubefile;
        cubefile.open(filename, std::fstream::out);

        const auto& steps = this->grid.steps();
        const auto& origin = this->grid.origin();
        const auto& step_sizes = this->grid.stepSizes();
        const auto& nuclei = molecule.nuclearFramework().nucleiAsVector();


        // Write the necessary header lines.

        // The first two lines are comment lines.
        cubefile << "COMMENT LINE -- GAUSSIAN Cube file" << std::endl;
        cubefile << "COMMENT LINE -- OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" << std::endl;

        cubefile << std::scientific;

        // The next line has the number of atoms and the origin of the volumetric data.
        cubefile << nuclei.size() << " " << origin(0) << " " << origin(1) << " " << origin(2) << std::endl;

        // The next three lines give the number of voxels along the respective axes.
        // We're choosing the x-, y- and z-axes, and since the number of steps is positive, the units are Bohr.
        cubefile << steps[0] << " " << step_sizes[0] << " " << 0.0 << " " << 0.0 << std::endl;
        cubefile << steps[1] << " " << 0.0 << " " << step_sizes[1] << " " << 0.0 << std::endl;
        cubefile << steps[2] << " " << 0.0 << " " << 0.0 << " " << step_sizes[2] << std::endl;
        for (const auto& nucleus : nuclei) {
            cubefile << nucleus.charge() << " " << 0.0 << " " << nucleus.position()(0) << " " << nucleus.position()(1) << " " << nucleus.position()(2) << std::endl;
        }


        // Write the values of the scalar function.
        size_t index = 0;
        this->grid.forEach([&index, &cubefile, this](const size_t i, const size_t j, const size_t k) {
            cubefile << this->values[index] << " ";  // write one value

            // There can only be 5 values on one line.
            if (k % 6 == 5) {
                cubefile << std::endl;
            }

            if (j == 0) {
                cubefile << std::endl;
            }

            index++;  // move to the next value
        });

        cubefile.close();
    }
};


}  // namespace GQCP
