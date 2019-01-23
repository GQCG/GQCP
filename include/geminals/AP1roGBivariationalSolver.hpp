#ifndef AP1roGBivariationalSolver_hpp
#define AP1roGBivariationalSolver_hpp


#include "geminals/BaseAP1roGSolver.hpp"


namespace GQCP {


/**
 *  A class that is able to solve the AP1roG bivariational equations
 */
class AP1roGBivariationalSolver : public BaseAP1roGSolver {
public:
    // CONSTRUCTORS
    /**
     *  @param N_P          the number of electrons
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *  @param G            the initial guess for the AP1roG gemial coefficients
     */
    AP1roGBivariationalSolver(size_t N_P, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G);

    /**
     *  @param N_P          the number of electrons
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGBivariationalSolver(size_t N_P, const HamiltonianParameters& ham_par);

    /**
     *  @param molecule     the molecule used for the AP1roG calculation
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *  @param G            the initial guess for the AP1roG gemial coefficients
     */
    AP1roGBivariationalSolver(const Molecule& molecule, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G);

    /**
     *  @param molecule     the molecule used for the AP1roG calculation
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGBivariationalSolver(const Molecule& molecule, const HamiltonianParameters& ham_par);


    // PUBLIC METHODS
    void solve() override;
};


}  // namespace GQCP


#endif /* AP1roGBivariationalSolver_hpp */
