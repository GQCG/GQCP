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
#include "Basis/BasisFunction.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param lc       a linear combination of CartesianGTOs
 */
BasisFunction::BasisFunction(const Base& lc) :
    Base(lc)
{
    this->N_total = this->calculateNormalizationFactor();
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the total normalization factor of this basis function
 */
double BasisFunction::calculateNormalizationFactor() const {
    return 0.0;
}


}  // namespace GQCP
