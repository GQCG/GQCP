#include "TwoElectronOperator.hpp"



#include <iostream>

namespace GQCG {


/**
 *  Constructor based on a given @param tensor
 */
TwoElectronOperator::TwoElectronOperator(const Eigen::Tensor<double, 4>& tensor) :
    tensor (tensor)
{
    // Check if the given tensor is 'square'
    auto dims = tensor.dimensions();
    if ((dims[0] != dims[1]) || (dims[1] != dims[2]) || (dims[2] != dims[3]) ) {
        throw std::invalid_argument("The given tensor should have equal dimensions in every rank.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Transform the two-electron integrals using the transformation matrix @param T
 *
 *  Note that the transformation matrix @param T is used as
 *      b' = b T ,
 *  in which the basis functions are collected as elements of a row vector b
 */
void TwoElectronOperator::transform(const Eigen::MatrixXd& T) {

}

/**
 *  Rotate the two-electron integrals using a unitary rotation matrix @param U
 *
 *  Note that the rotation matrix @param U is used as
 *      b' = b U ,
 *  in which the basis functions are collected as elements of a row vector b.
 */
void TwoElectronOperator::rotate(const Eigen::MatrixXd& U) {

}

/**
 *  Rotate the two-electron integrals using the unitary Jacobi rotation matrix U constructed from the @param jacobi_rotation_parameters
 *
 *  Note that
 *      - the rotation matrix @param U is used as
 *          b' = b U ,
 *        in which the basis functions are collected as elements of a row vector b.
 *      - we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
 */
void TwoElectronOperator::rotate(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters) {

}


}  // namespace GQCG
