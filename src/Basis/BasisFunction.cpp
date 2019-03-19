#include "Basis/BasisFunction.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param lc       a linear combination of CartesianGTOs
 */
BasisFunction::BasisFunction(const Base& lc) :
    Base(lc)
{
    this->N_total = this->calculateNormalizationFactor();
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the total normalization factor of this basis function
 */
double BasisFunction::calculateNormalizationFactor() const {
    return 0.0;
}


}  // namespace GQCP
