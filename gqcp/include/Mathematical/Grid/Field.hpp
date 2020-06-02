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


#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/miscellaneous.hpp"
#include "Utilities/type_traits.hpp"

#include <vector>


namespace GQCP {


/**
 *  A set of function values corresponding to points in space.
 * 
 *  @tparam T           the type of the evaluated function values
 */
template <typename T>
class Field {
private:
    std::vector<T> m_values;  // the evaluated function values, in the order of the grid's loop


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  The memberwise constructor.
     * 
     *  @param values           the evaluated function values, in the order of the grid's loop
     */
    Field(const std::vector<T>& values) :
        m_values {values} {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Parse a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/) for scalar field values. The grid-associated information is discarded.
     *
     *  @param filename                 the name of the cubefile
     * 
     *  @note This named constructor is only enabled for Field<double>.
     */
    template <typename Z = T>
    static enable_if_t<std::is_same<Z, double>::value, Field<double>> ReadCubeFile(const std::string& filename) {

        // Prepare the input file.
        std::ifstream input_file_stream = validateAndOpen(filename, "cube");


        // Skip the first two comment lines.
        std::string line;
        std::getline(input_file_stream, line);
        std::getline(input_file_stream, line);


        // Figure out the number of atoms, in order to ignore that information later on.
        std::getline(input_file_stream, line);

        // Split the line on any whitespace or tabs.
        std::vector<std::string> splitted_line;  // create a container for the line to be split in

        boost::trim_if(line, boost::is_any_of(" \t"));
        boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        const auto number_of_atoms = std::stoi(splitted_line[0]);


        // The next three lines contain the number of steps.
        size_t grid_size {1};  // will contain the total grid size after parsing
        for (size_t i = 0; i < 3; i++) {
            std::getline(input_file_stream, line);

            // Split the line on any whitespace or tabs.
            boost::trim_if(line, boost::is_any_of(" \t"));
            boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

            // The first column contains the number of steps.
            grid_size *= static_cast<size_t>(std::stoll(splitted_line[0]));
        }


        // Skip the lines containing atomic information.
        for (size_t i = 0; i < number_of_atoms; i++) {
            std::getline(input_file_stream, line);
        }


        // Finally, read in all field values.
        std::vector<double> field_values;
        field_values.reserve(grid_size);  // allocate memory for all the grid points

        while (std::getline(input_file_stream, line)) {

            // Split the line on any whitespace or tabs.
            boost::trim_if(line, boost::is_any_of(" \t"));
            boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

            // Add the field values from this line to the total array.
            for (const auto& value : splitted_line) {
                field_values.push_back(std::stod(value));  // no worries about push_back, since we have reserved enough memory up-front
            }
        }

        input_file_stream.close();

        return Field<double>(field_values);
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the size of this field, i.e. the number of field values
     */
    size_t size() const { return this->m_values.size(); }

    /**
     *  Access one of the field's values.
     * 
     *  @param index                the index of the function value
     * 
     *  @return a read-only field value, corresponding to the given index
     */
    const T& value(const size_t index) const { return this->m_values[index]; }

    /**
     *  Access one of the field's values.
     * 
     *  @param index                the index of the function value
     * 
     *  @return a writable field value, corresponding to the given index
     */
    T& value(const size_t index) { return this->m_values[index]; }

    /**
     *  @return the evaluated function values, in the order of the grid's loop
     */
    const std::vector<T>& values() const { return this->m_values; }
};


/*
 *  Convenience aliases for fields.
 */

template <typename Scalar>
using MatrixField = Field<Matrix<Scalar, 3, 3>>;

template <typename Scalar>
using VectorField = Field<Vector<Scalar, 3>>;


}  // namespace GQCP
