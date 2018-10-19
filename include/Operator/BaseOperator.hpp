#ifndef GQCP_BASEOPERATOR_HPP
#define GQCP_BASEOPERATOR_HPP


#include <Eigen/Dense>

#include "JacobiRotationParameters.hpp"



namespace GQCP {


/**
 *  A base class for the representation of operators in an orbital basis
 */
class BaseOperator {
protected:
    const size_t dim;  // dimension of the matrix representation of the operator


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param dimension
     */
    explicit BaseOperator(size_t dimension);


    // GETTERS
    size_t get_dim() { return this->dim; }


    // PUBLIC METHODS
    /**
     *  Transform the matrix representation of an operator using the transformation matrix @param T
     *
     *  Note that the transformation matrix @param T is used as
     *      b' = b T ,
     *  in which the basis functions are collected as elements of a row vector b
     */
    virtual void transform(const Eigen::MatrixXd& T) = 0;

    /**
     *  Rotate the matrix representation of an operator using a unitary rotation matrix @param U
     *
     *  Note that the rotation matrix @param U is used as
     *      b' = b U ,
     *  in which the basis functions are collected as elements of a row vector b.
     */
    virtual void rotate(const Eigen::MatrixXd& U) = 0;

    /**
     *  Rotate the matrix representation of an operator using the unitary Jacobi rotation matrix U constructed from the @param jacobi_rotation_parameters
     *
     *  Note that
     *      - the rotation matrix @param U is used as
     *          b' = b U ,
     *        in which the basis functions are collected as elements of a row vector b.
     *      - we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    virtual void rotate(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters) = 0;
};



}  // namespace GQCP


#endif  // GQCP_BASEOPERATOR_HPP
