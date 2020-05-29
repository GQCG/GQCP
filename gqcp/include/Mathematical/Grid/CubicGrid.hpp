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


#include "Mathematical/Grid/Field.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Mathematical/AbstractFunction/ScalarFunction.hpp"

#include <array>
#include <functional>
#include <type_traits>


namespace GQCP {


/**
 *  A grid type whose points are on a regular cubic lattice.
 */
class CubicGrid {
private:
    Vector<double, 3> m_origin;        // the origin of the grid
    std::array<size_t, 3> m_steps;     // the number of steps in the x, y, z-directions
    std::array<double, 3> step_sizes;  // the step sizes in the x, y, z-directions


public:
    // CONSTRUCTORS

    /**
     *  @param origin               the origin of the grid
     *  @param steps                the number of steps in the x, y, z-directions
     *  @param step_sizes           the step sizes in the x, y, z-directions
     */
    CubicGrid(const Vector<double, 3>& origin, const std::array<size_t, 3>& steps, const std::array<double, 3>& step_sizes);


    // PUBLIC METHODS

    /**
     *  Evaluate a scalar function on every point of this grid.
     * 
     *  @param scalar_function          the scalar functions whose values should be evaluated
     * 
     *  @return a field with the calculated evaluations
     */
    template <typename Valued>
    Field<Valued, CubicGrid> evaluate(const ScalarFunction<Valued, double, 3>& scalar_function) const {

        std::vector<Valued> values;  // the evaluated values of the scalar function
        values.reserve(this->numberOfPoints());

        this->forEach([&values, &scalar_function](const Vector<double, 3>& r) {
            values.push_back(scalar_function(r));
        });

        return Field<Valued, CubicGrid> {values, *this};
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
     *  @return the number of steps in the x, y, z-directions
     */
    const std::array<size_t, 3>& steps() const { return this->m_steps; }

    /**
     *  @return the step sizes in the x, y, z-directions
     */
    const std::array<double, 3>& stepSizes() const { return this->step_sizes; }
};


}  // namespace GQCP
