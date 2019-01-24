#ifndef BivariationalCoefficients_hpp
#define BivariationalCoefficients_hpp


#include "geminals/AP1roGVariables.hpp"


namespace GQCP {


/**
 *  A struct that holds the solutions (q0, q_i^a) to the bivariational equations
 */
struct BivariationalCoefficients {
    double q0;
    AP1roGVariables q;
};


}  // namespace GQCP



#endif /* BivariationalCoefficients_hpp */
