// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2019  the GQCG developers
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
#include "QCMethod/RMP2/RMP2.hpp"


namespace GQCP {


/**
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param molecule                     the molecule for which the energy correction should be calculated
 *  @param rhf_parameters               the converged solution to the RHF SCF equations
 *
 *  @return the RMP2 energy correction
 */
double calculateRMP2EnergyCorrection(const SQHamiltonian<double>& sq_hamiltonian, const Molecule& molecule, const QCModel::RHF<double>& rhf_parameters) {

    const size_t N = molecule.numberOfElectrons();
    const size_t K = sq_hamiltonian.dimension();

    const size_t HOMO_index = QCModel::RHF<double>::HOMOIndex(N);
    const size_t LUMO_index = QCModel::RHF<double>::LUMOIndex(K, N);

    const auto& g = sq_hamiltonian.twoElectron().parameters();

    double E = 0.0;
    //  loop over all occupied orbitals (0 <= HOMO )
    for (size_t i = 0; i <= HOMO_index; i++) {
        for (size_t j = 0; j <= HOMO_index; j++) {

            //  loop over all virtual orbitals (LUMO < K)
            for (size_t a = LUMO_index; a < K; a++) {
                for (size_t b = LUMO_index; b < K; b++) {

                    double epsilon_a = rhf_parameters.orbitalEnergy(a);
                    double epsilon_b = rhf_parameters.orbitalEnergy(b);
                    double epsilon_i = rhf_parameters.orbitalEnergy(i);
                    double epsilon_j = rhf_parameters.orbitalEnergy(j);

                    E -= g(a, i, b, j) * (2 * g(i, a, j, b) - g(i, b, j, a)) / (epsilon_a + epsilon_b - epsilon_i - epsilon_j);
                }
            }  // end of summation over virtual orbitals
        }
    }  // end of summation over occupied orbitals

    return E;
}


}  // namespace GQCP
