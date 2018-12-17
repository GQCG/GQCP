// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef GQCP_JACOBIROTATIONPARAMETERS_HPP
#define GQCP_JACOBIROTATIONPARAMETERS_HPP


#include <stdlib.h>
#include <ostream>

#include <Eigen/Dense>


namespace GQCP {


/**
 *  A class that holds the parameters that define a Jacobi rotation
 */
class JacobiRotationParameters {
private:
    size_t p;  // p > q, starts from 0
    size_t q;  // starts from 0
    double angle;  // expressed in radians

public:
    // CONSTRUCTORS
    /**
     *  @param p        the index of the first rotated orbital
     *  @param q        the index of the second rotated orbital
     *  @param angle    the angle of rotation, in radians
     */
    JacobiRotationParameters(size_t p, size_t q, double angle);


    // OPERATORS
    /**
     *  @param os                               the output stream which the parameters should be concatenated to
     *  @param jacobi_rotation_parameters       the parameters that should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const JacobiRotationParameters& jacobi_rotation_parameters);


    // GETTERS
    size_t get_p() const { return this->p; }
    size_t get_q() const { return this->q; }
    double get_angle() const { return this->angle; }
};



}  // namespace GQCP


#endif  // GQCP_JACOBIROTATIONPARAMETERS_HPP
