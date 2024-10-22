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


#include "Mathematical/Functions/Function.hpp"
#include "Mathematical/Grid/Field.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Molecule/Molecule.hpp"

#include <array>
#include <functional>
#include <numeric>
#include <type_traits>


namespace GQCP {


/**
 *  A grid type whose points are on a regular cubic lattice.
 */
class CubicGrid {
private:
    Vector<double, 3> m_origin;              // the origin of the grid
    std::array<size_t, 3> numbers_of_steps;  // the number of steps in the x, y, z-directions
    std::array<double, 3> step_sizes;        // the step sizes in the x, y, z-directions


public:
    // CONSTRUCTORS

    /**
     *  @param origin                       the origin of the grid
     *  @param numbers_of_steps             the number of steps in the x, y, z-directions
     *  @param step_sizes                   the step sizes in the x, y, z-directions
     */
    CubicGrid(const Vector<double, 3>& origin, const std::array<size_t, 3>& numbers_of_steps, const std::array<double, 3>& step_sizes);


    // NAMED CONSTRUCTORS

    /**
     *  Create a cubic ground that is centered around the given point. Equal numbers of steps and step sizes are taken in the x-, y- and z-directions.
     * 
     *  @param point                        the point around which the grid should be centered
     *  @param number_of_steps              the number of steps that should be taken in the x-, y- and z-directions
     *  @param step_size                    the size of the steps that should be taken in the x-, y- and z-directions
     */
    static CubicGrid Centered(const Vector<double, 3>& point, const size_t number_of_steps, const double step_size);

    /**
     *  Parse a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/). The values for the contained scalar field are ignored.
     *
     *  @param filename                 the name of the cubefile
     * 
     *  @note The Cube file is assumed to have grid axes oriented along the x-, y-, and z-axes.
     */
    static CubicGrid ReadCubeFile(const std::string& filename);

    /**
     *  Parse an .rgrid-file and create the CubicGrid that is contained in it. The values for the scalar field or vector field are ignored.
     * 
     *  @param filename             the name of the .igrid-file
     * 
     *  @note A regular grid (.rgrid) file is a headerless file and contains the following data:
     *      - Each row relates to one grid point, where the fastest changing values are z > y > x.
     *      - Column specification:
     *          - Column 1: The index from 1 to the number of grid points
     *          - Columns 2-4: The position of the grid point: x, y, and z
     *          - Optional: Column 5 or columns 5-7: 1 value for a scalar field, 3 values for a vector field
     */
    static CubicGrid ReadRegularGridFile(const std::string& filename);


    // PUBLIC METHODS

    /**
     *  Evaluate a function on every point of this grid.
     * 
     *  @param function          The function whose values should be evaluated.
     * 
     *  @return A field with the calculated evaluations.
     */
    template <typename OutputType>
    Field<OutputType> evaluate(const Function<OutputType, Vector<double, 3>>& function) const {

        std::vector<OutputType> values;
        values.reserve(this->numberOfPoints());

        this->forEach([&values, &function](const Vector<double, 3>& r) {
            values.push_back(function(r));
        });

        return Field<OutputType>(values);
    }


    /**
     *  Integrate a Field over this grid.
     * 
     *  @param field            the field that should be integrated, i.e. provided as the integrand
     * 
     *  @return the value of the integral
     */
    template <typename T>
    T integrate(const Field<T>& field) const {

        // Using field.value(0) as the initial value makes sure that the initialization is already correct: e.g. when a Vector is initialized, it doesn't automatically have the required size due to the default constructor
        auto result = std::accumulate(field.values().begin(), field.values().end(), field.value(0));
        result -= field.value(0);  // subtract the initial value of result

        // For cubic grids, the weight of each point is the voxel volume, so we'll have to multiply the final result by this.
        return result * this->voxelVolume();
    }


    /**
     *  Loop over the points of this grid by index number.
     * 
     *  @param callback         the function you would like to apply to each incoming (i,j,k)-tuple of numbers of steps taken in the x,y,z-direction.
     */
    void forEach(const std::function<void(const size_t, const size_t, const size_t)>& callback) const;

    /**
     *  Loop over the points of this grid by position (relative to the origin of this grid).
     * 
     *  @param callback         the function you would like to apply to each incoming position vector
     */
    void forEach(const std::function<void(const Vector<double, 3>&)>& callback) const;

    /**
     *  @return the number of points that are in this grid
     */
    size_t numberOfPoints() const;

    /**
     *  @return the origin of this grid
     */
    const Vector<double, 3>& origin() const { return this->m_origin; }

    /**
     *  @param i        the number of steps taken in the x-direction
     *  @param j        the number of steps taken in the y-direction
     *  @param k        the number of steps taken in the z-direction
     *
     *  @return the position vector associated to the given indices
     */
    Vector<double, 3> position(const size_t i, const size_t j, const size_t k) const;

    /**
     *  @return a vector of the points that are described by this grid
     */
    std::vector<Vector<double, 3>> points() const;

    /**
     *  @param axis         0, 1, 2 representing the x-, y-, or z-axis
     * 
     *  @return the number of steps that can be taken in the direction of the specified axis
     */
    size_t numbersOfSteps(const size_t axis) const { return this->numbers_of_steps[axis]; }

    /**
     *  @return the number of steps in the x, y, z-directions
     */
    const std::array<size_t, 3>& numbersOfSteps() const { return this->numbers_of_steps; }

    /**
     *  @param axis         0, 1, 2 representing the x-, y-, or z-axis
     * 
     *  @return the step size that is taken in the direction of the specified axis
     */
    double stepSize(const size_t axis) const { return this->step_sizes[axis]; }

    /**
     *  @return the step sizes in the x, y, z-directions
     */
    const std::array<double, 3>& stepSizes() const { return this->step_sizes; }

    /**
     *  @return the total volume that is contained in this grid
     */
    double totalVolume() const { return this->voxelVolume() * this->numberOfPoints(); }

    /**
     *  Write a field's values to a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/).
     *
     *  @param scalar_field             the scalar field that should be written to the cubefile
     *  @param filename                 the name of the cubefile that has to be generated
     *  @param molecule                 the molecule that should be placed in the cubefile
     */
    void writeToCubeFile(const Field<double>& scalar_field, const std::string& filename, const Molecule& molecule) const;

    /**
     *  @return the volume of one voxel in this grid
     */
    double voxelVolume() const;
};


}  // namespace GQCP
