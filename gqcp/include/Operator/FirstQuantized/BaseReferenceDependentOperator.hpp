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


namespace GQCP {


/**
 *  A base class/interface used to represent first-quantized operators that are dependent on a point of reference.
 * 
 *  @example Some examples of reference-dependent operators are the multipole operators and the angular momentum operator
 */
class BaseReferenceDependentOperator {
private:
    // The point that is used as a reference to define the operator.
    Vector<double, 3> m_reference;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a `BaseReferenceDependentOperator` using a point of reference.
     * 
     *  @param reference            The point that is used as a reference to define the operator.
     */
    BaseReferenceDependentOperator(const Vector<double, 3>& reference = Vector<double, 3>::Zero());


    /*
     *  MARK: Destructor
     */

    // Make the destructor pure virtual in order to make this class abstract.
    virtual ~BaseReferenceDependentOperator() = 0;


    /*
     *  MARK: Reference
     */

    /**
     *  @return The point that is used as a reference to define the operator.
     */
    const Vector<double, 3>& reference() const { return this->m_reference; }
};


}  // namespace GQCP
