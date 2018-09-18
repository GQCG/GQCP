#ifndef GQCG_ONEELECTRONOPERATOR_HPP
#define GQCG_ONEELECTRONOPERATOR_HPP


#include <Eigen/Dense>

#include "BaseOperator.hpp"


namespace GQCG {



class OneElectronOperator : public BaseOperator {
private:
    Eigen::MatrixXd matrix;  // the matrix representation of the one-electron operator


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param matrix
     */
    explicit OneElectronOperator(const Eigen::MatrixXd& matrix);


    // GETTERS
    Eigen::MatrixXd get_matrix_representation() const { return this->matrix; }


    // PUBLIC METHODS
    /**
     *  Transform the matrix representation of a one-electron operator using the transformation matrix @param T
     *
     *  Note that the transformation matrix @param T is used as
     *      b' = b T ,
     *  in which the basis functions are collected as elements of a row vector b
     */
    void transform(const Eigen::MatrixXd& T) override;

    /**
     *  Rotate the matrix representation of a one-electron operator using a unitary rotation matrix @param U
     *
     *  Note that the rotation matrix @param U is used as
     *      b' = b U ,
     *  in which the basis functions are collected as elements of a row vector b.
     */
    void rotate(const Eigen::MatrixXd& U) override;

    /**
     *  Rotate the matrix representation of a one-electron operator using the unitary Jacobi rotation matrix U constructed from the @param jacobi_rotation_parameters
     *
     *  Note that
     *      - the rotation matrix @param U is used as
     *          b' = b U ,
     *        in which the basis functions are collected as elements of a row vector b.
     *      - we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters) override;
};



}  // namespace GQCG



#endif  // GQCG_ONEELECTRONOPERATOR_HPP
