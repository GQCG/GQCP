#ifndef GQCP_CONSTRAINT_HPP
#define GQCP_CONSTRAINT_HPP

#include "OneElectronOperator.hpp"
#include "TwoElectronOperator.hpp"

namespace GQCP {



struct Constraint {

    OneElectronOperator* one_op = nullptr;
    TwoElectronOperator* two_op = nullptr;

};


}

#endif  // GQCP_CONSTRAINT_HPP
