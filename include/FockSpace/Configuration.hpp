#ifndef GQCP_CONFIGURATION_HPP
#define GQCP_CONFIGURATION_HPP


#include "ONV.hpp"


namespace GQCP {


/**
 *  Struct containing an electron distribution
 *  represented by the alpha and beta ONV.
 */
struct Configuration {
    ONV onv_alpha;
    ONV onv_beta;
};


}  // namespace GQCP


#endif  // GQCP_CONFIGURATION_HPP
