#include "OneElectronOperator.hpp"

#include <stdexcept>


namespace GQCG {


/*
 *  CONSTRUCTOR
 */

/**
 *  Constructor based on a given @param matrix
 */
OneElectronOperator::OneElectronOperator(const Eigen::MatrixXd& matrix) :
    matrix (matrix)
{
    // Check if the one-electron integrals are represented as a square matrix
    if (matrix.cols() != matrix.rows()) {
        throw std::invalid_argument("One-electron integrals have to be represented as a square matrix.");
    }
}


/*
 *  PUBLIC METHODS
 */

/**
 *  Transform the one-electron integrals using the transformation matrix @param T
 *
 *  Note that the transformation matrix @param T is used as
 *      b' = b T ,
 *  in which the basis functions are collected as elements of a row vector b
 */
void OneElectronOperator::transform(const Eigen::MatrixXd& T) {

}


/**
 *  Rotate the one-electron integrals using a unitary rotation matrix @param U
 *
 *  Note that the rotation matrix @param U is used as
 *      b' = b U ,
 *  in which the basis functions are collected as elements of a row vector b.
 */
void OneElectronOperator::rotate(const Eigen::MatrixXd& U) {

}


/**
 *  Rotate the one-electron integrals using the unitary Jacobi rotation matrix U constructed from the @param jacobi_rotation_parameters
 *
 *  Note that
 *      - the rotation matrix @param U is used as
 *          b' = b U ,
 *        in which the basis functions are collected as elements of a row vector b.
 *      - we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
 */
void OneElectronOperator::rotate(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters) {

}



}  // namespace GQCG
