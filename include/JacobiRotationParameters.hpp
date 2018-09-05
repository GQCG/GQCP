#ifndef GQCG_JACOBIROTATIONPARAMETERS_HPP
#define GQCG_JACOBIROTATIONPARAMETERS_HPP


#include <stdlib.h>


namespace GQCG {


/**
 *  A class that holds @member p, @member q and @member angle to define a Jacobi rotation
 *
 *  Note that:
 *      - @param p must always be smaller than @param q
 *      - @member angle is expressed in radians
 */
class JacobiRotationParameters {
private:
    const size_t p;  // p < q
    const size_t q;
    const double angle;

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
};



}  // namespace GQCG


#endif  // GQCG_JACOBIROTATIONPARAMETERS_HPP
