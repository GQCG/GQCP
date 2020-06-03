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
    // CONSTRUCTORS

    /**
     *  A memberwise constructor.
     * 
     *  @param points               the grid points
     *  @param weights              a 1-D array containing the weights for each of the grid points
     */
    WeightedGrid(const std::vector<Vector<double, 3>>& points, const ArrayX<double>& weights);


    // NAMED CONSTRUCTORS

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
    static WeightedGrid ReadIntegrationGridFile(const std::string& filename);


    // PUBLIC METHODS

    /**
     *  Integrate a Field over this grid.
     * 
     *  @param field            the field that should be integrated, i.e. provided as the integrand
     * 
     *  @return the value of the integral
     */
    template <typename T>
    T integrate(const Field<T>& field) const {

        // A KISS-implementation of integration: multiplying every field value by the weight of the associated grid point.
        auto result = field.value(0) * this->weight(0);  // this makes sure that 'result' is already initialized to the correct type, e.g. when a Vector is initialized, it doesn't automatically have the required size
        for (size_t i = 1; i < this->size(); i++) {
            result += field.value(i) * this->weight(i);
        }

        return result;
    }


    /**
     *  @return the number of grid points/weights
     */
    size_t numberOfPoints() const { return this->m_points.size(); }

    /**
     *  Access one of the grid's points.
     * 
     *  @param index                the index of the grid point
     * 
     *  @return a read-only grid point, corresponding to the given index
     */
    const Vector<double, 3>& point(const size_t index) const { return this->m_points[index]; }

    /**
     *  @return the grid points
     */
    const std::vector<Vector<double, 3>>& points() const { return this->m_points; }

    /**
     *  @return the size of the grid, i.e. the number of grid points/weights
     */
    size_t size() const { return this->numberOfPoints(); }

    /**
     *  Access one of the grid's weights.
     * 
     *  @param index                the index of the weight
     * 
     *  @return a read-only grid weight, corresponding to the given index
     */
    double weight(const size_t index) const { return this->m_weights(index); }

    /**
     *  @return a 1-D array containing the weights for each of the grid points
     */
    const ArrayX<double>& weights() const { return this->m_weights; };
};


}  // namespace GQCP