#include "OrbitalOptimization/AP1roGLagrangianNewtonOrbitalOptimizer.hpp"

#include "Geminals/AP1roGLagrangianOptimizer.hpp"
#include "Geminals/AP1roG.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param N_P          the number of electron pairs
 *  @param G            the initial guess for the AP1roG gemial coefficients
 */
AP1roGLagrangianNewtonOrbitalOptimizer::AP1roGLagrangianNewtonOrbitalOptimizer(size_t N_P, const AP1roGGeminalCoefficients& G, std::shared_ptr<NewtonOrbitalOptimizationOptions> oo_options) :
    N_P (N_P),
    G (G),
    QCMethodNewtonOrbitalOptimizer(oo_options)
{}


/**
 *  @param N_P          the number of electron pairs
 *  @param K            the number of spatial orbitals
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGLagrangianNewtonOrbitalOptimizer::AP1roGLagrangianNewtonOrbitalOptimizer(size_t N_P, size_t K, std::shared_ptr<NewtonOrbitalOptimizationOptions> oo_options) : 
    AP1roGLagrangianNewtonOrbitalOptimizer(N_P, AP1roGGeminalCoefficients(N_P, K), std::move(oo_options))
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer
 * 
 *  In the case of this uncoupled AP1roG Lagrangian orbital optimizer, the PSEs are re-solved in every iteration using the current orbitals
 */
void AP1roGLagrangianNewtonOrbitalOptimizer::prepareDMCalculation(const HamiltonianParameters<double>& ham_par) {

    // Solve the AP1roG PSEs and determine the Lagrangian multipliers
    AP1roGLagrangianOptimizer lagrangian_optimizer (this->N_P, ham_par, this->G);
    lagrangian_optimizer.solve();
    this->E = lagrangian_optimizer.get_electronic_energy();
    this->G = lagrangian_optimizer.get_geminal_coefficients();
    this->multipliers = lagrangian_optimizer.get_multipliers();

    std::cout << "Current energy: " << this->E << std::endl << std::endl;
    std::cout << "Current amplitues: " << std::endl << this->G.asMatrix() << std::endl << std::endl;
}


/**
 *  @return the current 1-DM
 */
OneRDM<double> AP1roGLagrangianNewtonOrbitalOptimizer::calculate1RDM() const {
    return GQCP::calculate1RDM(this->G, this->multipliers);
}


/**
 *  @return the current 2-DM
 */
TwoRDM<double> AP1roGLagrangianNewtonOrbitalOptimizer::calculate2RDM() const {
    return GQCP::calculate2RDM(this->G, this->multipliers);
}


/**
 *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
 * 
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return the new full set orbital generators, including the redundant parameters
 */
OrbitalRotationGenerators AP1roGLagrangianNewtonOrbitalOptimizer::calculateNewFullOrbitalGenerators(const HamiltonianParameters<double>& ham_par) const {
    std::cout << "Current kappa: " << std::endl << this->calculateNewFreeOrbitalGenerators(ham_par).asVector() << std::endl << std::endl;
    return this->calculateNewFreeOrbitalGenerators(ham_par);  // no extra step necessary
}


}  // namespace GQCP
