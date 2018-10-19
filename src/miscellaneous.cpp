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
