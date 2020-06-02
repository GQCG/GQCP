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

#include "Mathematical/Grid/CubicGrid.hpp"

#include "Utilities/miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param origin               the origin of the grid
 *  @param steps                the number of steps in the x, y, z-directions
 *  @param step_sizes           the step sizes in the x, y, z-directions
 */
CubicGrid::CubicGrid(const Vector<double, 3>& origin, const std::array<size_t, 3>& steps, const std::array<double, 3>& step_sizes) :
    m_origin {origin},
    number_of_steps {steps},
    step_sizes {step_sizes} {}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Parse a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/). The values for the contained scalar field are ignored.
 *
 *  @param filename                 the name of the cubefile
 * 
 *  @note The Cube file is assumed to have grid axes oriented along the x-, y-, and z-axes.
 */
CubicGrid CubicGrid::ReadCubeFile(const std::string& filename) {

    // Prepare the input file stream and some variables.
    std::ifstream input_file_stream = validateAndOpen(filename, "cube");

    Vector<double, 3> origin = Vector<double, 3>::Zero();
    std::array<double, 3> step_sizes {0.0, 0.0, 0.0};
    std::array<size_t, 3> number_of_steps {0, 0, 0};


    // Do the actual parsing.
    std::string line;


    // Skip the first two comment lines.
    std::getline(input_file_stream, line);
    std::getline(input_file_stream, line);


    // Read in the origin of the grid.
    std::getline(input_file_stream, line);

    // Split the line on any whitespace or tabs.
    std::vector<std::string> splitted_line;  // create a container for the line to be split in

    boost::trim_if(line, boost::is_any_of(" \t"));
    boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

    const auto x = std::stod(splitted_line[1]);
    const auto y = std::stod(splitted_line[2]);
    const auto z = std::stod(splitted_line[3]);

    origin << x, y, z;


    // The next three lines contain the step sizes and the number of steps.
    for (size_t i = 0; i < 3; i++) {
        std::getline(input_file_stream, line);

        // Split the line on any whitespace or tabs.
        boost::trim_if(line, boost::is_any_of(" \t"));
        boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        // The first column contains the number of steps.
        number_of_steps[i] = static_cast<size_t>(std::stoll(splitted_line[0]));

        // The next three columns contain the step size.
        step_sizes[i] = std::stod(splitted_line[i + 1]);
    }

    input_file_stream.close();

    return CubicGrid(origin, number_of_steps, step_sizes);
}


/**
 *  Parse an .rgrid-file and create the CubicGrid that is contained in it. The values for the scalar field or vector field are ignored.
 * 
 *  @param filename             the name of the .igrid-file
 * 
 *  @note An integration grid (.igrid) file is a headerless file and contains the following data:
 *      - Each row relates to one grid point, where the fastest changing values are z > y > x.
 *      - Column specification:
 *          - Column 1: The index from 1 to the number of grid points
 *          - Columns 2-4: The position of the grid point: x, y, and z
 *          - Optional: Column 5 or columns 5-7: 1 value for a scalar field, 3 values for a vector field
 */
CubicGrid CubicGrid::ReadRegularGridFile(const std::string& filename) {

    // Prepare the input file stream and some variables.
    std::ifstream input_file_stream = validateAndOpen(filename, "rgrid");

    Vector<double, 3> origin = Vector<double, 3>::Zero();
    std::array<double, 3> step_sizes {0.0, 0.0, 0.0};
    std::array<size_t, 3> number_of_steps {0, 0, 0};


    // Do the actual parsing.
    std::string line;


    // We'll treat the first line as the origin.
    std::getline(input_file_stream, line);

    // Split the line on any whitespace or tabs.
    std::vector<std::string> splitted_line;  // create a container for the line to be split in

    boost::trim_if(line, boost::is_any_of(" \t"));
    boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

    // Read the coordinates of the grid point.
    auto x = std::stod(splitted_line[1]);
    auto y = std::stod(splitted_line[2]);
    auto z = std::stod(splitted_line[3]);

    origin << x, y, z;


    // Continue parsing, by figuring out the step sizes and the number of steps in each Cartesian direction.
    // Assume that the fastest varying axis are z > y > x.

    // If we read one line, we can figure out the step size in the z-direction.
    std::getline(input_file_stream, line);
    boost::trim_if(line, boost::is_any_of(" \t"));
    boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

    z = std::stod(splitted_line[3]);

    step_sizes[2] = z - origin(2);


    // Keep reading lines until the y-coordinate changes.
    while (std::getline(input_file_stream, line)) {

        // Split the line on any whitespace or tabs.
        boost::trim_if(line, boost::is_any_of(" \t"));
        boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        // Read the y-coordinate of the grid point. If it has changed, count the number of steps taken in the z-dimension
        y = std::stod(splitted_line[2]);
        if (y != origin(1)) {
            // Read the index column and fill in the number of steps in the z-direction.
            const auto index = static_cast<size_t>(std::stoll(splitted_line[0]));

            number_of_steps[2] = index - 1;

            // Since the y-coordinate changed, we can figure out the step size in the y-direction.
            step_sizes[1] = y - origin(1);
            break;
        }
    }


    // Keep reading lines until the x-coordinate changes.
    while (std::getline(input_file_stream, line)) {

        // Split the line on any whitespace or tabs.
        boost::trim_if(line, boost::is_any_of(" \t"));
        boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        // Read the x-coordinate of the grid point. If it has changed, count the number of steps taken in the y-direction.
        x = std::stod(splitted_line[1]);
        if (x != origin(0)) {
            // Read the index column and fill in the number of steps in the y-direction.
            const auto index = static_cast<size_t>(std::stoll(splitted_line[0]));

            number_of_steps[1] = (index - 1) / number_of_steps[2];

            // Since the x-coordinate changed, we can figure out the step size in the x-direction.
            step_sizes[0] = x - origin(0);
            break;
        }
    }


    // Read until the end of the file to figure out the number of steps taken in the x-direction.
    size_t final_index;  // will eventually contain the final index
    while (std::getline(input_file_stream, line)) {

        // Split the line on any whitespace or tabs.
        boost::trim_if(line, boost::is_any_of(" \t"));
        boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        final_index = static_cast<size_t>(std::stoll(splitted_line[0]));  // will eventually contain the final index
    }
    number_of_steps[0] = final_index / (number_of_steps[1] * number_of_steps[2]);


    // We're done parsing now.
    input_file_stream.close();

    return CubicGrid(origin, number_of_steps, step_sizes);
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of points that are in this grid
 */
size_t CubicGrid::numberOfPoints() const {

    return (this->step_sizes[0] + 1) * (this->step_sizes[1] + 1) * (this->step_sizes[2] + 1);
}


/**
 *  Loop over the points of this grid by index number.
 * 
 *  @param callback         the function you would like to apply to each incoming (i,j,k)-tuple of numbers of steps taken in the x,y,z-direction.
 */
void CubicGrid::forEach(const std::function<void(const size_t, const size_t, const size_t)>& callback) const {

    for (size_t i = 0; i < this->number_of_steps[0]; i++) {
        for (size_t j = 0; j < this->number_of_steps[1]; j++) {
            for (size_t k = 0; k < this->number_of_steps[2]; k++) {
                callback(i, j, k);
            }
        }
    }
}


/**
 *  Loop over the points of this grid by position (relative to the origin of this grid).
 * 
 *  @param callback         the function you would like to apply to each incoming position vector
 */
void CubicGrid::forEach(const std::function<void(const Vector<double, 3>&)>& callback) const {

    const auto this_copy = *this;
    this->forEach([this_copy, callback](const size_t i, const size_t j, const size_t k) {
        const auto position = this_copy.position(i, j, k);
        callback(position);
    });
}


/**
 *  @param i        the number of steps taken in the x-direction
 *  @param j        the number of steps taken in the y-direction
 *  @param k        the number of steps taken in the z-direction
 *
 *  @return the position vector associated to the given indices
 */
Vector<double, 3> CubicGrid::position(const size_t i, const size_t j, const size_t k) const {

    const double x = this->m_origin(0) + i * this->step_sizes[0];
    const double y = this->m_origin(1) + j * this->step_sizes[1];
    const double z = this->m_origin(2) + k * this->step_sizes[2];

    return Vector<double, 3>(x, y, z);
}


/**
 *  Write a field's values to a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/).
 *
 *  @param scalar_field             the scalar field that should be written to the cubefile
 *  @param filename                 the name of the cubefile that has to be generated
 *  @param molecule                 the molecule that should be placed in the cubefile
 */
void CubicGrid::writeToCubeFile(const Field<double>& scalar_field, const std::string& filename, const Molecule& molecule) const {

    // Prepare some variables.
    std::ofstream cubefile;
    cubefile.open(filename, std::fstream::out);

    const auto& steps = this->numberOfSteps();
    const auto& origin = this->origin();
    const auto& step_sizes = this->stepSizes();
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
    this->forEach([&index, &cubefile, &scalar_field](const size_t i, const size_t j, const size_t k) {
        cubefile << scalar_field.value(index) << " ";  // write one value

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


/**
 *  @return the volume of the voxels in this grid
 */
double CubicGrid::voxelVolume() const {

    double volume = 1.0;
    for (const auto& step_size : this->step_sizes) {
        volume *= step_size;
    }
    return volume;
}


}  // namespace GQCP
