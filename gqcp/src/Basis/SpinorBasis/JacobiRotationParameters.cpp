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

#include "Basis/SpinorBasis/JacobiRotationParameters.hpp"

#include <stdexcept>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param p        the index of the first rotated orbital
 *  @param q        the index of the second rotated orbital
 *  @param angle    the angle of rotation, in radians
 */
JacobiRotationParameters::JacobiRotationParameters(const size_t p, const size_t q, const double angle) :
    m_p {p},
    m_q {q},
    m_angle {angle} {

    // Check if p > q
    if (this->m_p <= this->m_q) {
        throw std::invalid_argument("JacobiRotationParameters::JacobiRotationParameters(size_t, size_t, double): Can't construct a JacobiRotationParameter with p < q.");
    }
}


/**
 *  Default constructor
 */
JacobiRotationParameters::JacobiRotationParameters() :
    JacobiRotationParameters(1, 0, 0.0) {}


/*
 *  OPERATORS
 */
/**
 *  @param os                               the output stream which the jacobi rotation parameters should be concatenated to
 *  @param jacobi_rotation_parameters       the parameters that should be concatenated to the output stream
 *
 *  @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const JacobiRotationParameters& jacobi_rotation_parameters) {

    os << "p: " << jacobi_rotation_parameters.p() << ", q: " << jacobi_rotation_parameters.q() << ", angle: " << jacobi_rotation_parameters.angle();
    return os;
}


}  // namespace GQCP
