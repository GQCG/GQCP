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


#include "Utilities/Eigen.hpp"

#include <cstdlib>
#include <ostream>


namespace GQCP {


/**
 *  A type that encapsulates the parameters that define a real Jacobi rotation.
 */
class JacobiRotationParameters {
private:
    size_t m_p;      // p > q, starts from 0
    size_t m_q;      // starts from 0
    double m_angle;  // expressed in radians

public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param p        the index of the first rotated orbital
     *  @param q        the index of the second rotated orbital
     *  @param angle    the angle of rotation, in radians
     */
    JacobiRotationParameters(const size_t p, const size_t q, const double angle);

    /**
     *  Default constructor setting everything to zero
     */
    JacobiRotationParameters();


    /*
     *  MARK: Operators
     */

    /**
     *  @param os                               the output stream which the parameters should be concatenated to
     *  @param jacobi_rotation_parameters       the parameters that should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const JacobiRotationParameters& jacobi_rotation_parameters);


    /*
     *  MARK: General information
     */

    /**
     *  @return the angle of rotation, in radians
     */
    double angle() const { return this->m_angle; }

    /**
     *  @return the index of the first rotated orbital
     */
    size_t p() const { return this->m_p; }

    /**
     *  @return the index of the second rotated orbital
     */
    size_t q() const { return this->m_q; }


    /*
     *  MARK: Conversions
     */

    /**
     *  @return This as an Eigen::JacobiRotation.
     */
    Eigen::JacobiRotation<double> Eigen() const { return Eigen::JacobiRotation<double> {std::cos(this->angle()), std::sin(this->angle())}; }
};


}  // namespace GQCP
