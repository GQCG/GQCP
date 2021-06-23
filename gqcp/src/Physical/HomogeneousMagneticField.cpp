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

#include "Physical/HomogeneousMagneticField.hpp"


namespace GQCP {

/*
 *  MARK: Constructors
 */

/**
 *  Initialize a `HomogeneousMagneticField` from a field strength and gauge origin.
 * 
 *  @param B            The field strength.
 *  @param G            The gauge origin.
 */
HomogeneousMagneticField::HomogeneousMagneticField(const Vector<double, 3>& B, const Vector<double, 3>& G) :
    B {B},
    G {G} {}


/*
 *  MARK: Physics
 */

/**
 *  Calculate the vector potential at a point in space.
 * 
 *  @param r        A point in space.
 * 
 *  @return The vector potential evaluated at the given point.
 */
Vector<double, 3> HomogeneousMagneticField::vectorPotentialAt(const Vector<double, 3>& r) const {

    return 0.5 * this->strength().cross(r - this->gaugeOrigin());
}


}  // namespace GQCP
