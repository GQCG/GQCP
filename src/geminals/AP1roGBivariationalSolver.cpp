#include "geminals/AP1roGBivariationalSolver.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 */
AP1roGBivariationalSolver::AP1roGBivariationalSolver(size_t N_P, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G) :
    BaseAP1roGSolver(N_P, ham_par, G)
{}


/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGBivariationalSolver::AP1roGBivariationalSolver(size_t N_P, const HamiltonianParameters& ham_par) :
    BaseAP1roGSolver(N_P, ham_par)
{}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 */
AP1roGBivariationalSolver::AP1roGBivariationalSolver(const Molecule& molecule, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G) :
    BaseAP1roGSolver(molecule, ham_par, G)
{}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGBivariationalSolver::AP1roGBivariationalSolver(const Molecule& molecule, const HamiltonianParameters& ham_par) :
    BaseAP1roGSolver(molecule, ham_par)
{}



/*
 *  PUBLIC METHODS
 */
void AP1roGBivariationalSolver::solve() {
    // not implemented yet
}


}  // namespace GQCP
