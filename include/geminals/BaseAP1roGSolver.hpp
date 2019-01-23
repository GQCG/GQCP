#ifndef BaseAP1roGSolver_hpp
#define BaseAP1roGSolver_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "geminals/AP1roGGeminalCoefficients.hpp"


namespace GQCP {


/**
 *  A base class for solvers using the AP1roG wave function
 */
class BaseAP1roGSolver {
protected:
    size_t K;  // the number of spatial orbitals
    size_t N_P;  // the number of electron pairs
    double electronic_energy;  // the converged electronic energy

    AP1roGGeminalCoefficients geminal_coefficients;  // the converged geminal coefficients

    HamiltonianParameters ham_par;


public:
    // CONSTRUCTORS
    /**
     *  @param N_P          the number of electrons
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *  @param G            the initial guess for the AP1roG gemial coefficients
     */
    BaseAP1roGSolver(size_t N_P, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G);

    /**
     *  @param N_P          the number of electrons
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *
     *  The initial guess for the geminal coefficients is zero
     */
    BaseAP1roGSolver(size_t N_P, const HamiltonianParameters& ham_par);

    /**
     *  @param molecule     the molecule used for the AP1roG calculation
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *  @param G            the initial guess for the AP1roG gemial coefficients
     */
    BaseAP1roGSolver(const Molecule& molecule, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G);

    /**
     *  @param molecule     the molecule used for the AP1roG calculation
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *
     *  The initial guess for the geminal coefficients is zero
     */
    BaseAP1roGSolver(const Molecule& molecule, const HamiltonianParameters& ham_par);


    // DESTRUCTOR
    virtual ~BaseAP1roGSolver();


    // GETTERS
    double get_electronic_energy() const { return this->electronic_energy; }
    const AP1roGGeminalCoefficients& get_geminal_coefficients() const { return this->geminal_coefficients; }
    const HamiltonianParameters& get_hamiltonian_parameters() const { return this->ham_par; }


    // PUBLIC METHODS
    /**
     *  The actual 'solving' step
     */
    virtual void solve() = 0;
};


}  // namespace GQCP


#endif /* BaseAP1roGSolver_hpp */
