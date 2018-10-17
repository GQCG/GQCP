#include "AP1roG/AP1roG.hpp"


namespace GQCG {

/*
 *  CONSTRUCTORS
 */
/**
 *  Default constructor setting everything to zero
 */
AP1roG::AP1roG() :
    geminal_coefficients (AP1roGGeminalCoefficients()),
    electronic_energy (0.0)
{}


/**
 *  Constructor based on given @param geminal_coefficients and @param electronic_energy
 */
AP1roG::AP1roG(const GQCG::AP1roGGeminalCoefficients& geminal_coefficients, double electronic_energy) :
    geminal_coefficients (geminal_coefficients),
    electronic_energy (electronic_energy)
{}



/*
 *  HELPER FUNCTIONS
 */
/**
 *  Calculate the AP1roG energy given AP1roG geminal coefficients @param G and Hamiltonian parameters @param ham_par
 */
double calculateAP1roGEnergy(const GQCG::AP1roGGeminalCoefficients& G, const GQCG::HamiltonianParameters& ham_par) {

    Eigen::MatrixXd h_SO = ham_par.get_h().get_matrix_representation();
    Eigen::Tensor<double, 4> g_SO = ham_par.get_g().get_matrix_representation();


    // KISS implementation of the AP1roG energy
    double E = 0.0;
    for (size_t j = 0; j < G.get_N_P(); j++) {
        E += 2 * h_SO(j,j);

        for (size_t k = 0; k < G.get_N_P(); k++) {
            E += 2 * g_SO(k,k,j,j) - g_SO(k,j,j,k);
        }

        for (size_t b = G.get_N_P(); b < G.get_K(); b++) {
            E += g_SO(j,b,j,b) * G(j,b);
        }
    }

    return E;
}


}  // namespace GQCG
