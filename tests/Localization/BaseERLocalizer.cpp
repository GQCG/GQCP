#include "Localization/BaseERLocalizer.hpp"


namespace GQCP {



/*
 *  CONSTRUCTORS
 */
/**
 *  @param N_P                              the number of electron pairs
 *  @param threshold                        the threshold for maximization on subsequent localization indices
 *  @param maximum_number_of_iterations     the maximum number of iterations for the localization algorithm
 */
BaseERLocalizer::BaseERLocalizer(size_t N_P, double threshold, size_t maximum_number_of_iterations) :
    N_P (N_P),
    threshold (threshold),
    maximum_number_of_iterations (maximum_number_of_iterations)
{}



}  // namespace GQCP
