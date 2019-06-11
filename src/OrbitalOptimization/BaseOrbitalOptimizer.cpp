#include "OrbitalOptimization/BaseOrbitalOptimizer.hpp"


namespace GQCP {


/* 
 *  CONSTRUCTORS
 */

/**
 *  @param oo_options       the orbital optimization options that should be used for the orbital optimization algorithm
 */
BaseOrbitalOptimizer::BaseOrbitalOptimizer(const OrbitalOptimizationOptions& oo_options) :
    oo_options (oo_options)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Optimize the Hamiltonian parameters by subsequently
 *      - checking for convergence (see checkForConvergence())
 *      - rotating the Hamiltonian parameters with a newly found rotation matrix (see calculateNewRotationMatrix())
 * 
 *  @param ham_par      the initial (guess for the) Hamiltonian parameters
 */
void BaseOrbitalOptimizer::optimize(HamiltonianParameters<double>& ham_par) {

    if (!ham_par.areOrbitalsOrthonormal()) {
        throw std::invalid_argument("BaseOrbitalOptimizer::optimize(HamiltonianParameters<double>&): The given Hamiltonian parameters do not belong to an orthonormal basis.");
    }

    size_t number_of_oo_iterations {0};
    while (this->prepareConvergenceChecking(ham_par), !this->checkForConvergence(ham_par)) {  // result of the comma operator is the second operand, if not converged
        this->prepareRotationMatrixCalculation(ham_par);
        auto U = this->calculateNewRotationMatrix(ham_par);
        ham_par.rotate(U);

        number_of_oo_iterations++;
        if (number_of_oo_iterations > this->oo_options.maximum_number_of_iterations) {
            throw std::invalid_argument("BaseOrbitalOptimizer::optimize(HamiltonianParameters<double>&): The orbital optimization procedure did not converge in the given amount of iterations.");
        }
    }

    this->is_converged = true;
}


}  // namespace GQCP