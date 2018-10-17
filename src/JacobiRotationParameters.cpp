#include "JacobiRotationParameters.hpp"

#include <stdexcept>


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param p, @param q and a @param angle expressed in radians
 */
JacobiRotationParameters::JacobiRotationParameters(size_t p, size_t q, double angle) :
    p (p),
    q (q),
    angle (angle)
{
    // Check if p > q
    if (this->p <= this->q) {
        throw std::invalid_argument("Can't construct a JacobiRotationParameter with p < q.");
    }
}


/*
 *  OPERATORS
 */
/**
 *  Overloading of operator<< for GQCG::JacobiRotationParameters to be used with streams
 */
std::ostream& operator<<(std::ostream& os, const GQCG::JacobiRotationParameters& jacobi_rotation_parameters) {

    os << "p: " << jacobi_rotation_parameters.p << ", q: " << jacobi_rotation_parameters.q << ", angle: " << jacobi_rotation_parameters.angle;
    return os;
}



}  // namespace GQCG
