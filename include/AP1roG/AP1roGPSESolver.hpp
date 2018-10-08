#ifndef AP1roGPSESolver_hpp
#define AP1roGPSESolver_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "AP1roG/AP1roGGeminalCoefficients.hpp"
#include "AP1roG/AP1roG.hpp"

namespace GQCG {

/**
 *
 */
class AP1roGPSESolver {
private:
    const size_t K;  // the number of special orbitals
    const size_t N_P;  // the number of electron pairs
    const GQCG::AP1roGGeminalCoefficients initial_geminal_coefficients;
    
    GQCG::HamiltonianParameters ham_par;

    GQCG::AP1roG solution;


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param molecule and Hamiltonian parameters @param ham_par
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGPSESolver(const GQCG::Molecule& molecule, const GQCG::HamiltonianParameters& ham_par);
};


}  // namespace GQCG



#endif /* AP1roGPSESolver_hpp */
