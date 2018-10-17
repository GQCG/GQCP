#ifndef GQCG_ONEELECTRONOPERATOR_HPP
#define GQCG_ONEELECTRONOPERATOR_HPP


#include <Eigen/Dense>

#include "BaseOperator.hpp"


namespace GQCG {



/**
 *  A class that holds the matrix representation of a one-electron operator in an orbital basis
 */
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
    double get(size_t p, size_t q) { return this->matrix(p,q); }

    
    // OPERATORS
    /**
     *  @return the sum of two OneElectronOperators, i.e. a OneElectronOperator whose matrix representation is the sum
     *  of the two matrix representations of the given OneElectronOperators
     */
    GQCG::OneElectronOperator operator+(const GQCG::OneElectronOperator& other);


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


    // FRIEND CLASSES
    friend class HamiltonianParameters;
};



}  // namespace GQCG



#endif  // GQCG_ONEELECTRONOPERATOR_HPP
