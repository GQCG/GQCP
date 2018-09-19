#ifndef GQCG_MISCELLANEOUS_HPP
#define GQCG_MISCELLANEOUS_HPP


#include <stdlib.h>

#include <Eigen/Dense>

#include "JacobiRotationParameters.hpp"


namespace GQCG {


/**
 *  @return the @param M-dimensional Jacobi rotation matrix given a set of JacobiRotationParameters
 *
 *  Note that:
 *      - we work with the (cos, sin, -sin, cos) definition of the Jacobi rotation matrix
 */
Eigen::MatrixXd jacobiRotationMatrix(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters, size_t M);



}  // namespace GQCG



#endif  // GQCG_MISCELLANEOUS_HPP
