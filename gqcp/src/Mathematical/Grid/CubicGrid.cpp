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
    m_steps {steps},
    step_sizes {step_sizes} {}


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
void CubicGrid::loop(const std::function<void(const size_t, const size_t, const size_t)>& callback) const {

    for (size_t i = 0; i < this->m_steps[0]; i++) {
        for (size_t j = 0; j < this->m_steps[1]; j++) {
            for (size_t k = 0; k < this->m_steps[2]; k++) {
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
void CubicGrid::loop(const std::function<void(const Vector<double, 3>&)>& callback) const {

    const auto this_copy = *this;
    this->loop([this_copy, callback](const size_t i, const size_t j, const size_t k) {
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


}  // namespace GQCP
