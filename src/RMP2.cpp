#include "RMP2.hpp"



namespace GQCG {


/**
 *  Calculate and @return the RMP2 energy correction
 */
double calculateRMP2EnergyCorrection(const GQCG::HamiltonianParameters& ham_par, const GQCG::Molecule& molecule) {

    Eigen::Tensor<double, 4> g = ham_par.g.get_matrix_representation();
    auto K = g.dimension(0);  // TODO: change to ham_par.dim


    size_t LUMO_index = this->rhf.LUMOIndex();
    size_t HOMO_index = this->rhf.HOMOIndex();

    double E = 0.0;
    //  loop over all occupied orbitals (0 <= HOMO )
    for (size_t i = 0; i <= HOMO_index; i++) {
        for (size_t j = 0; j <= HOMO_index; j++) {

            //  loop over all virtual orbitals (LUMO < K)
            for (size_t a = LUMO_index; a < K; a++) {
                for (size_t b = LUMO_index; b < K; b++) {

                    double epsilon_a = this->rhf.get_orbital_energies(a);
                    double epsilon_b = this->rhf.get_orbital_energies(b);
                    double epsilon_i = this->rhf.get_orbital_energies(i);
                    double epsilon_j = this->rhf.get_orbital_energies(j);

                    E -= g(a,i,b,j) * (2 * g(i,a,j,b) - g(i,b,j,a)) / (epsilon_a + epsilon_b - epsilon_i - epsilon_j);
                }
            }  // end of summation over virtual orbitals


        }
    }  // end of summation over occupied orbitals

    return E;
}

}  // namespace GQCG
