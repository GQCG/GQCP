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
#pragma once

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
    m_origin (origin),
    m_steps (steps),
    step_sizes (step_sizes)
{}



/*
 *  PUBLIC METHODS
 */

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
