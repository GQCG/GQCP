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
#ifndef OrbitalOptimizationOptions_hpp
#define OrbitalOptimizationOptions_hpp


#include <cstdlib>


namespace GQCP {


/**
 *  A struct that holds options for orbital optimization
 */
struct OrbitalOptimizationOptions {
    double convergence_threshold = 1.0e-08;
    size_t maximum_number_of_iterations = 128;
};


}  // namespace GQCP


#endif /* OrbitalOptimizationOptions_hpp */
