#ifndef AP1roGBivariationalSolver_hpp
#define AP1roGBivariationalSolver_hpp


#include "geminals/BaseAP1roGSolver.hpp"
#include "geminals/AP1roGVariables.hpp"


namespace GQCP {


/**
 *  A struct that holds the solutions (q0, q_i^a) to the bivariational equations
 */
struct BivariationalCoefficients {
    double q0;
    AP1roGVariables q;
};



/**
 *  A class that is able to solve the AP1roG bivariational equations
 */
class AP1roGBivariationalSolver : public BaseAP1roGSolver {
private:
    BivariationalCoefficients bivariational_coefficients;  // the determined bivariational coefficients


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


    // GETTERS
    const BivariationalCoefficients& get_bivariational_coefficients() const { return this->bivariational_coefficients; }


    // PUBLIC METHODS
    void solve() override;
};


}  // namespace GQCP


#endif /* AP1roGBivariationalSolver_hpp */
