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
#include "miscellaneous.hpp"


namespace GQCP {

/**
 *  @return the @param M-dimensional Jacobi rotation matrix given a set of JacobiRotationParameters
 *
 *  Note that:
 *      - we work with the (cos, sin, -sin, cos) definition of the Jacobi rotation matrix
 */
Eigen::MatrixXd jacobiRotationMatrix(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters, size_t M) {

    double c = std::cos(jacobi_rotation_parameters.angle);
    double s = std::sin(jacobi_rotation_parameters.angle);

    // We'll start the construction with an identity matrix
    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(M, M);

    // And apply the Jacobi rotation as J = I * jacobi_rotation (cfr. B' = B T)
    J.applyOnTheRight(jacobi_rotation_parameters.p, jacobi_rotation_parameters.q, Eigen::JacobiRotation<double> (c, s));
    return J;
}


}  // namespace GQCP
