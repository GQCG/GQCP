#include "Mathematical/Optimization/UnalteringHessianModifier.hpp"


namespace GQCP {


/**
 *  @param hessian      the current indefinite Hessian
 * 
 *  @return the given Hessian, i.e. do not alter the current hessian
 */
SquareMatrix<double> UnalteringHessianModifier::operator()(const SquareMatrix<double>& hessian) {
    return hessian;
}


}  // namespace GQCP
