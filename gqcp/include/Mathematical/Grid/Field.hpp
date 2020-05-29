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
     *  PUBLIC METHODS
     */

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
