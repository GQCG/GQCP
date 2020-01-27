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
#pragma once


#include "ONV.hpp"


namespace GQCP {


/**
 *  A struct that represents the combination of an alpha and beta ONV
 */
struct Configuration {
    ONV onv_alpha;
    ONV onv_beta;

    /**
     *  @param delimiter        the delimiter between bits
     * 
     *  @return a compact representation of the configuration as string:
     *   e.g: 1100 | 1010 => 2110
     */
    std::string spinSummedRepresentation(const std::string& delimiter = "") const;
};


}  // namespace GQCP
