#include "AP1roGPSESolver.hpp"


namespace GQCG {


AP1roGPSESolver::AP1roGPSESolver(const GQCG::Molecule& molecule, const GQCG::HamiltonianParameters& ham_par) :
    K (ham_par.K),
    ham_par (ham_par),
    N_P (molecule.N / 2),
    initial_geminal_coefficients (GQCG::AP1roGGeminalCoefficients(this->N_P, this->K))
{
    // Check if we have an even number of electrons
    if ((molecule.N % 2) != 0) {
        throw std::invalid_argument("The given molecule has an odd number of electrons.");
    }
}


}  // namespace GQCG
