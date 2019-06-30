#include "Operator/FirstQuantized/BaseNuclearOperator.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param atoms                the atoms that represent the nuclear framework
 */
BaseNuclearOperator::BaseNuclearOperator(const std::vector<Atom>& atoms) :
    atoms (atoms)
{}


}  // namespace GQCP
