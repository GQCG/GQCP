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

#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"


namespace GQCP {


/*
 *  MARK: Operator value
 */

/**
 *  @return The scalar value of this nuclear repulsion operator.
 */
double NuclearRepulsionOperator::value() const {

    const auto& nuclei = this->nuclear_framework.nucleiAsVector();

    // Sum over every unique nucleus pair.
    double value {0.0};
    const auto n_nuclei = this->nuclearFramework().numberOfNuclei();
    for (size_t i = 0; i < n_nuclei; i++) {
        for (size_t j = i + 1; j < n_nuclei; j++) {
            const auto nucleus1 = nuclei[i];
            const auto nucleus2 = nuclei[j];

            // The internuclear repulsion energy (Coulomb) for every nucleus pair is Z1 * Z2 / |R1 - R2|.
            value += nucleus1.charge() * nucleus2.charge() / nucleus1.calculateDistanceWith(nucleus2);
        }
    }

    return value;
}


}  // namespace GQCP
