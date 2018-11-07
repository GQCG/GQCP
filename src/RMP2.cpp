// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "RMP2.hpp"


namespace GQCP {


/**
 *  @param ham_par      Hamiltonian parameters in an orthornomal orbital basis
 *  @param molecule     the molecule for which the energy correction should be calculated
 *  @param rhf          the converged solution to the RHF SCF equations
 *
 *  @return the RMP2 energy correction
 */
double calculateRMP2EnergyCorrection(const GQCP::HamiltonianParameters& ham_par, const GQCP::Molecule& molecule, const GQCP::RHF& rhf) {

    size_t N = molecule.get_N();
    size_t K = ham_par.get_K();

    size_t HOMO_index = GQCP::RHFHOMOIndex(N);
    size_t LUMO_index = GQCP::RHFLUMOIndex(K, N);


    GQCP::TwoElectronOperator g = ham_par.get_g();

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
