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


#include "Mathematical/Representation/Array.hpp"
#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/**
 *  A collection of points in 3D-space, with each point associated to a weight.
 */
class WeightedGrid {
private:
    ArrayX<double> m_weights;  // a 1-D array containing the weights for each of the grid points

    std::vector<Vector<double, 3>> m_points;  // the grid points


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  A memberwise constructor.
     * 
     *  @param points               the grid points
     *  @param weights              a 1-D array containing the weights for each of the grid points
     */
    WeightedGrid(const std::vector<Vector<double, 3>>& points, const ArrayX<double>& weights);


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the grid points
     */
    const std::vector<Vector<double, 3>>& points() const { return this->m_points; }

    /**
     *  @return a 1-D array containing the weights for each of the grid points
     */
    const ArrayX<double> weights() const { return this->m_weights; };
};


}  // namespace GQCP