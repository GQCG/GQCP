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


}  // namespace GQCP
