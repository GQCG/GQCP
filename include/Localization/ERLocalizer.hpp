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
#ifndef ERLocalizer_hpp
#define ERLocalizer_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCP {


/**
 *  A class that localizes a set of orthonormal orbitals according to the Edmiston-Ruedenberg localization index
 */
class ERLocalizer {
private:
    const size_t N_P;     // the number of electron pairs


public:
    // CONSTRUCTORS
    /**
     *  @param N_P        the number of electron pairs
     */
    ERLocalizer(size_t N_P);


    // PUBLIC METHODS
    /**
     *  @param ham_par      the Hamiltonian parameters that contain the two-electron integrals upon which the Edmiston-Ruedenberg localization index is calculated
     *
     *  @return the Edmiston-Ruedenberg localization index
     */
    double calculateLocalizationIndex(const GQCP::HamiltonianParameters& ham_par) const;
};


}  // namespace GQCP


#endif /* ERLocalizer_hpp */
