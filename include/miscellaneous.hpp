#ifndef GQCP_MISCELLANEOUS_HPP
#define GQCP_MISCELLANEOUS_HPP


#include <stdlib.h>

#include <Eigen/Dense>

#include "JacobiRotationParameters.hpp"


namespace GQCP {


/**
 *  @return the @param M-dimensional Jacobi rotation matrix given a set of JacobiRotationParameters
 *
 *  Note that:
 *      - we work with the (cos, sin, -sin, cos) definition of the Jacobi rotation matrix
 */
Eigen::MatrixXd jacobiRotationMatrix(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters, size_t M);



}  // namespace GQCP



#endif  // GQCP_MISCELLANEOUS_HPP
