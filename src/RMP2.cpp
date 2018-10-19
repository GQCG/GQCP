#include "RMP2.hpp"


namespace GQCP {


/**
 *  @return the RMP2 energy correction based on given @param Hamiltonian parameters ham_par, a given @param molecule and a converged solution @param rhf to the RHF SCF equations
 */
double calculateRMP2EnergyCorrection(const GQCP::HamiltonianParameters& ham_par, const GQCP::Molecule& molecule, const GQCP::RHF& rhf) {

    size_t N = molecule.get_N();

    Eigen::Tensor<double, 4> g = ham_par.get_g().get_matrix_representation();
    size_t K = g.dimension(0);  // TODO: change to ham_par.dim after rebase


    size_t HOMO_index = GQCP::RHFHOMOIndex(N);
    size_t LUMO_index = GQCP::RHFLUMOIndex(K, N);

    double E = 0.0;
    //  loop over all occupied orbitals (0 <= HOMO )
    for (size_t i = 0; i <= HOMO_index; i++) {
        for (size_t j = 0; j <= HOMO_index; j++) {

            //  loop over all virtual orbitals (LUMO < K)
            for (size_t a = LUMO_index; a < K; a++) {
                for (size_t b = LUMO_index; b < K; b++) {

                    double epsilon_a = rhf.get_orbital_energies(a);
                    double epsilon_b = rhf.get_orbital_energies(b);
                    double epsilon_i = rhf.get_orbital_energies(i);
                    double epsilon_j = rhf.get_orbital_energies(j);

                    E -= g(a,i,b,j) * (2 * g(i,a,j,b) - g(i,b,j,a)) / (epsilon_a + epsilon_b - epsilon_i - epsilon_j);
                }
            }  // end of summation over virtual orbitals

        }
    }  // end of summation over occupied orbitals

    return E;
}


}  // namespace GQCP
