#ifndef BaseERLocalizer_hpp
#define BaseERLocalizer_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCP {


/**
 *  A base class to implement the Edmiston-Ruedenberg localization method
 */
class BaseERLocalizer {
protected:
    const size_t N_P;  // the number of electron pairs
    const double threshold;  // the threshold for maximization on subsequent localization indices
    const size_t maximum_number_of_iterations;  // the maximum number of iterations for the localization algorithm

    bool is_converged = false;
    size_t iterations = 0;  // the number of iterations

public:
    // CONSTRUCTORS
    /**
     *  @param N_P                              the number of electron pairs
     *  @param threshold                        the threshold for maximization on subsequent localization indices
     *  @param maximum_number_of_iterations     the maximum number of iterations for the localization algorithm
     */
    BaseERLocalizer(size_t N_P, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);


    // PUBLIC METHODS
    /**
     *  Localize the Hamiltonian parameters by maximizing the Edmiston-Ruedenberg localization index
     *
     *  @param ham_par      the Hamiltonian parameters that should be localized
     */
    virtual void localize(GQCP::HamiltonianParameters& ham_par) = 0;
};



}  // namespace GQCP


#endif /* BaseERLocalizer_hpp */
