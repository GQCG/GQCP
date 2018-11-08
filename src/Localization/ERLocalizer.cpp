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
#include "Localization/ERLocalizer.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param N_P        the number of electron pairs
 */
ERLocalizer::ERLocalizer(size_t N_P) :
    N_P (N_P)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param ham_par      the Hamiltonian parameters that contain the two-electron integrals upon which the Edmiston-Ruedenberg localization index is calculated
 *
 *  @return the Edmiston-Ruedenberg localization index
 */
double ERLocalizer::calculateLocalizationIndex(const GQCP::HamiltonianParameters& ham_par) const {

    double localization_index = 0.0;
    auto g = ham_par.get_g();  // the two-electron integrals


    // TODO: when Eigen releases TensorTrace, use it here
    for (size_t i = 0; i < this->N_P; i++) {
        localization_index += g(i,i,i,i);
    }

    return localization_index;
}




}  // namespace GQCP
