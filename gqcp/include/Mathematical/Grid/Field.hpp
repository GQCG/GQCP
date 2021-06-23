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

#include <algorithm>
#include <functional>
#include <vector>


namespace GQCP {


/**
 *  A set of function values corresponding to points in space.
 * 
 *  @tparam T_          The type of the evaluated function values.
 */
template <typename T_>
class Field {
public:
    // The type of the evaluated function values.
    using T = T_;

private:
    // The evaluated function values, in the order of the grid's loop.
    std::vector<T> m_values;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  The memberwise constructor.
     * 
     *  @param values           The evaluated function values, in the order of the grid's loop.
     */
    Field(const std::vector<T>& values) :
        m_values {values} {}


    /*
     *  MARK: Named constructors.
     */

    /**
     *  Parse a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/) for scalar field values. The grid-associated information is discarded.
     *
     *  @param filename                 The name of the cubefile.
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


    /**
     *  Parse an `.igrid`- or `.rgrid`-file and create the field that is contained in it. The values for the grid are discarded.
     * 
     *  @param filename             The name of the `.igrid`-file.
     * 
     *  @note An integration grid (`.igrid`) file is a headerless file and contains the following data:
     *      - Each row relates to one grid point.
     *      - Column specification:
     *          - Column 1: The index from 1 to the number of grid points
     *          - Columns 2-4: The position of the grid point: x, y, and z
     *          - Optional: Column 5 or columns 5-7: 1 value for a scalar field, 3 values for a vector field
     *          - Column 5, 6 or 8: The integration weight associated to the grid point
     * 
     *  @note A regular grid (`.rgrid`) file is a headerless file and contains the following data:
     *      - Each row relates to one grid point, where the fastest changing values are z > y > x.
     *      - Column specification:
     *          - Column 1: The index from 1 to the number of grid points
     *          - Columns 2-4: The position of the grid point: x, y, and z
     *          - Optional: Column 5 or columns 5-7: 1 value for a scalar field, 3 values for a vector field
     */
    template <int N>
    static Field<Vector<double, N>> ReadGridFile(const std::string& filename) {

        // Find the extension of the given path (https://stackoverflow.com/a/51992).
        std::string filename_extension;  // will contain the extension of the given filename
        std::string::size_type idx = filename.rfind('.');

        if (idx != std::string::npos) {
            filename_extension = filename.substr(idx + 1);
        } else {
            throw std::invalid_argument("Field::ReadGridFile(const std::string&): I did not find an extension in your given file name.");
        }


        // Check if the user supplied an .rgrid- or .igrid-file.
        if ((filename_extension != "rgrid") && (filename_extension != "igrid")) {
            throw std::invalid_argument("Field::ReadGridFile(const std::string&): The given file is not an .igrid- or .rgrid-file.");
        }


        // If the filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
        std::ifstream input_file_stream {filename};
        if (!input_file_stream.good()) {
            throw std::invalid_argument("validateAndOpen(const std::string&, const std::string&): The provided file name is illegible. Maybe you specified a wrong path?");
        }


        // Do the actual parsing.
        std::vector<Vector<double, N>> field_values;

        std::string line;
        while (std::getline(input_file_stream, line)) {
            std::vector<std::string> splitted_line;  // create a container for the line to be split in

            // Split the line on any whitespace or tabs.
            boost::trim_if(line, boost::is_any_of(" \t"));
            boost::split(splitted_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

            // We can skip the first 4 columns (index and grid point coordinates).

            // Read the field value by converting the string(s) to a double and using Eigen::Map to instantiate the correct Vector.
            std::vector<double> field_value_vector;
            for (size_t index = 4; index < 4 + N; index++) {
                field_value_vector.push_back(std::stod(splitted_line[index]));
            }
            Eigen::Map<Eigen::VectorXd> field_value {field_value_vector.data(), N};
            field_values.push_back(Vector<double, N>(field_value));


            // We can skip the weight column as well.
        }
        input_file_stream.close();


        // Convert the std::vector of weights to an Array.
        return Field<Vector<double, N>>(field_values);
    }


    /*
     *  MARK: Canonical mathematical operators (https://en.cppreference.com/w/cpp/language/operators)
     */

    /**
     *  Add a field to this one by a point-wise addition.
     * 
     *  @param rhs          The right-hand side of the addition.
     * 
     *  @return A reference to the modified `this`.
     * 
     *  @note The fields are supposed to be defined on the same grid.
     */
    Field<T>& operator+=(const Field<T>& rhs) {

        std::transform(this->m_values.begin(), this->m_values.end(), rhs.values().begin(),
                       this->m_values.begin(), std::plus<T>());

        return *this;
    }


    /**
     *  Add two fields by a point-wise addition.
     * 
     *  @param lhs          The left-hand side of the addition.
     *  @param rhs          The right-hand side of the addition.
     * 
     *  @return The sum of the two fields.
     * 
     *  @note The fields are supposed to be defined on the same grid.
     */
    friend Field<T> operator+(Field<T> lhs, const Field<T>& rhs) {

        lhs += rhs;
        return lhs;
    }


    /**
     *  Negate the values of this field.
     * 
     *  @return A negated version of this field.
     */
    Field<T> operator-() const {

        std::vector<T> result;
        result.reserve(this->size());

        std::transform(this->m_values.begin(), this->m_values.end(),
                       std::back_inserter(result), std::negate<T>());

        return result;
    }


    /**
     *  Subtract another field from this field by point-wise subtraction.
     * 
     *  @param rhs              The right-hand side of the subtraction.
     * 
     *  @return A reference to the modified `this`.
     * 
     *  @note The fields are supposed to be defined on the same grid.
     */
    Field<T>& operator-=(const Field<T>& rhs) {

        // Define subtraction as addition with the negation.
        *this += (-rhs);
        return *this;
    }


    /**
     *  Subtract one field from another by point-wise subtraction.
     * 
     *  @param lhs          The left-hand side of the subtraction.
     *  @param rhs          The right-hand side of the subtraction.
     * 
     *  @return The difference of the two fields.
     * 
     *  @note The fields are supposed to be defined on the same grid.
     */
    friend Field<T> operator-(Field<T> lhs, const Field<T>& rhs) {

        lhs -= rhs;
        return lhs;
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return The size of this field, i.e. the number of field values.
     */
    size_t size() const { return this->m_values.size(); }


    /*
     *  MARK: Mappings
     */

    /**
     *  Apply a given function on each of this field's values, in-place.
     * 
     *  @param function             The function to be used on each of the values.
     */
    void map(const std::function<T(const T&)>& function) {

        std::transform(this->m_values.begin(), this->m_values.end(), this->m_values.begin(),
                       function);
    }


    /**
     *  Apply a given function on each of this field's values.
     * 
     *  @param function             The function to be used on each of the values.
     * 
     *  @return A field with the mapped values.
     */
    Field<T> mapped(const std::function<T(const T&)>& function) const {

        auto this_copy = *this;
        this_copy.map(function);

        return this_copy;
    }


    /*
     *  MARK: Access
     */

    /**
     *  Access one of the field's values.
     * 
     *  @param index                The index of the function value.
     * 
     *  @return A read-only field value, corresponding to the given index.
     */
    const T& value(const size_t index) const { return this->m_values[index]; }

    /**
     *  Access one of the field's values.
     * 
     *  @param index                The index of the function value.
     * 
     *  @return A writable field value, corresponding to the given index.
     */
    T& value(const size_t index) { return this->m_values[index]; }

    /**
     *  @return The evaluated function values, in the order of the grid's loop.
     */
    const std::vector<T>& values() const { return this->m_values; }
};


/*
 *  Convenience aliases for fields.
 */

template <typename Scalar>
using VectorField = Field<Vector<Scalar, 3>>;

template <typename Scalar>
using MatrixField = Field<Matrix<Scalar, 3, 3>>;


}  // namespace GQCP
