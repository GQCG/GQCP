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
#ifndef BivariationalCoefficients_hpp
#define BivariationalCoefficients_hpp


#include "geminals/AP1roGVariables.hpp"


namespace GQCP {


/**
 *  A struct that holds the solutions (q0, q_i^a) to the bivariational equations
 */
struct BivariationalCoefficients {
    double q0;
    AP1roGVariables q;
};


}  // namespace GQCP



#endif /* BivariationalCoefficients_hpp */
