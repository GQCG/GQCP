#ifndef GQCG_JACOBIROTATIONPARAMETERS_HPP
#define GQCG_JACOBIROTATIONPARAMETERS_HPP


#include <stdlib.h>

#include <Eigen/Dense>


namespace GQCG {


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


    // GETTERS
    size_t get_p() const { return this->p; }
    size_t get_q() const { return this->q; }
    double get_angle() const { return this->angle; }


    // FRIEND CLASSES
    friend class OneElectronOperator;
    friend class AP1roGJacobiOrbitalOptimizer;

    // FRIEND FUNCTIONS
    friend Eigen::MatrixXd jacobiRotationMatrix(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters, size_t M);
};



}  // namespace GQCG


#endif  // GQCG_JACOBIROTATIONPARAMETERS_HPP
