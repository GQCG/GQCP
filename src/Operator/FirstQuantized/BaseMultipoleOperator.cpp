#include "Operator/FirstQuantized/BaseMultipoleOperator.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param o        the origin of the multipole
 */
BaseMultipoleOperator::BaseMultipoleOperator(const Vector<double, 3>& o) :
    o (o)
{}


}  // namespace GQCP
