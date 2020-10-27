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

#include "Basis/Transformations/JacobiRotation.hpp"

#include <stdexcept>


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param p            The index of the first rotated orbital. (p > q)
 *  @param q            The index of the second rotated orbital. (p > q)
 *  @param angle        The angle of rotation, in radians.
 */
JacobiRotation::JacobiRotation(const size_t p, const size_t q, const double angle) :
    m_p {p},
    m_q {q},
    m_angle {angle} {

    // Check if p > q
    if (this->m_p <= this->m_q) {
        throw std::invalid_argument("JacobiRotation::JacobiRotation(size_t, size_t, double): Can't construct a JacobiRotationParameter with p < q.");
    }
}


/**
 *  The default constructor.
 */
JacobiRotation::JacobiRotation() :
    JacobiRotation(1, 0, 0.0) {}


/*
 *  MARK: Operators
 */

/**
 *  @param os                               The output stream which the parameters should be concatenated to.
 *  @param jacobi_rotation                  The Jacobi rotation.
 *
 *  @return A reference to the updated output stream.
 */
std::ostream& operator<<(std::ostream& os, const JacobiRotation& jacobi_rotation) {

    os << "p: " << jacobi_rotation.p() << ", q: " << jacobi_rotation.q() << ", angle: " << jacobi_rotation.angle();
    return os;
}


}  // namespace GQCP
