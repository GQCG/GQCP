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

#include "Mathematical/Grid/WeightedGrid.hpp"

#include "Utilities/miscellaneous.hpp"

#include <boost/algorithm/string.hpp>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  A memberwise constructor.
 * 
 *  @param points               the grid points
 *  @param weights              a 1-D array containing the weights for each of the grid points
 */
WeightedGrid::WeightedGrid(const std::vector<Vector<double, 3>>& points, const ArrayX<double>& weights) :
    m_weights {weights},
    m_points {points} {

    if (this->m_weights.size() != this->m_points.size()) {
        throw std::invalid_argument("WeightedGrid(const std::vector<Vector<double, 3>>&, const ArrayX<double>&: The number of weights does not match the number of points.");
    }
}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Convert a cubic grid into a weighted grid, where the weights are all equal to the CubicGrid's voxel volume.
 * 
 *  @param cubic_grid               the cubic grid
 */
WeightedGrid WeightedGrid::FromCubicGrid(const CubicGrid& cubic_grid) {

    // Generate the grid points and associated weights.
    const auto points = cubic_grid.points();
    ArrayX<double> weights = ArrayX<double>::Constant(cubic_grid.numberOfPoints(), cubic_grid.voxelVolume());

    return WeightedGrid(points, weights);
}


/**
 *  Parse an .igrid-file and create the WeightedGrid that is contained in it. The values for the scalar field or vector field are discarded.
 * 
 *  @param filename             the name of the .igrid-file
 * 
 *  @note An integration grid (.igrid) file is a headerless file and contains the following data:
 *      - Each row relates to one grid point.
 *      - Column specification:
 *          - Column 1: The index from 1 to the number of grid points
 *          - Columns 2-4: The position of the grid point: x, y, and z
 *          - Optional: Column 5 or columns 5-7: 1 value for a scalar field, 3 values for a vector field
 *          - Column 5, 6 or 8: The integration weight associated to the grid point
 */
WeightedGrid WeightedGrid::ReadIntegrationGridFile(const std::string& filename) {

    // Prepare the input file and the containers for the grid points and associated weights.
    std::ifstream input_file_stream = validateAndOpen(filename, "igrid");

    std::vector<Vector<double, 3>> points;
    std::vector<double> weights;


    // Do the actual parsing.
    std::string line;
    while (std::getline(input_file_stream, line)) {
        std::vector<std::string> splitted_line;  // create a container for the line to be split in

        // Split the line on any whitespace or tabs.
        boost::trim_if(line, boost::is_any_of(" \t"));
        boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        // We can skip the index column.

        // Read the coordinates of the grid point.
        const auto x = std::stod(splitted_line[1]);
        const auto y = std::stod(splitted_line[2]);
        const auto z = std::stod(splitted_line[3]);

        points.emplace_back(x, y, z);


        // Read the associated weight.
        // Since we're discarding the info of the scalar or vector field contained in the file, we can be safe by retrieving the weight as the last column instead of specifying a hard-coded position.
        const auto weight = std::stod(splitted_line.back());

        weights.push_back(weight);
    }
    input_file_stream.close();


    // Convert the std::vector of weights to an Array.
    Eigen::Map<Eigen::ArrayXd> weights_array {weights.data(), static_cast<long>(weights.size())};  // need static_cast to Eigen::Index
    return WeightedGrid(points, weights_array);
}


}  // namespace GQCP
