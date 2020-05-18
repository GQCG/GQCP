// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "QCMethod/RMP2/RMP2.hpp"


namespace GQCP {


/**
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param rhf_parameters               the converged solution to the RHF SCF equations
 *
 *  @return the RMP2 energy correction
 */
double calculateRMP2EnergyCorrection(const SQHamiltonian<double>& sq_hamiltonian, const QCModel::RHF<double>& rhf_parameters) {

    // Prepare some variables.
    const auto orbital_space = rhf_parameters.orbitalSpace();
    const auto& g = sq_hamiltonian.twoElectron().parameters();


    double E = 0.0;
    for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
        double epsilon_i = rhf_parameters.orbitalEnergy(i);

        for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
            double epsilon_j = rhf_parameters.orbitalEnergy(j);

            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                double epsilon_a = rhf_parameters.orbitalEnergy(a);

                for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                    double epsilon_b = rhf_parameters.orbitalEnergy(b);


                    E -= g(a, i, b, j) * (2 * g(i, a, j, b) - g(i, b, j, a)) / (epsilon_a + epsilon_b - epsilon_i - epsilon_j);
                }
            }  // end of summation over virtual orbitals
        }
    }  // end of summation over occupied orbitals

    return E;
}


}  // namespace GQCP
