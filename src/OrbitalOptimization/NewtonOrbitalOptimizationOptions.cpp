#include "OrbitalOptimization/NewtonOrbitalOptimizationOptions.hpp"


namespace GQCP {


/**
 *  @param hessian_modifier         the modifier functor that should be used when an indefinite Hessian is encountered
 */
NewtonOrbitalOptimizationOptions::NewtonOrbitalOptimizationOptions(std::shared_ptr<BaseHessianModifier> hessian_modifier, const double convergence_threshold, const double maximum_number_of_iterations) :
    hessian_modifier(std::move(hessian_modifier)),
    OrbitalOptimizationOptions(convergence_threshold, maximum_number_of_iterations)
{}


}  // namespace GQCP
