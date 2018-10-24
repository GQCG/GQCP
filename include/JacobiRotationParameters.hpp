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
 *  A class that holds @member p, @member q and @member angle to define a Jacobi rotation
 *
 *  Note that:
 *      - @member p and @member q are indices that start from 0
 *      - @member p must always be larger than @member q
 *      - @member angle is expressed in radians
 */
class JacobiRotationParameters {
private:
    size_t p;  // p > q
    size_t q;
    double angle;

public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param p, @param q and a @param angle expressed in radians
     */
    JacobiRotationParameters(size_t p, size_t q, double angle);


    // OPERATORS
    /**
     *  Overloading of operator<< for GQCP::JacobiRotationParameters to be used with streams
     */
    friend std::ostream& operator<<(std::ostream& os, const GQCP::JacobiRotationParameters& jacobi_rotation_parameters);


    // GETTERS
    size_t get_p() const { return this->p; }
    size_t get_q() const { return this->q; }
    double get_angle() const { return this->angle; }
};



}  // namespace GQCP


#endif  // GQCP_JACOBIROTATIONPARAMETERS_HPP
