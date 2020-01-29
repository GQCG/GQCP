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


#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Optimization/Eigenproblem/DenseDiagonalization.hpp"


namespace GQCP {
namespace EigenproblemSolver {


/**
 *  @return an algorithm that can diagonalize a dense matrix
 */
Algorithm<EigenproblemEnvironment> Dense(const size_t number_of_requested_eigenpairs = 1) {

    // Our dense eigenproblem solver is just a wrapper around Eigen's routines.
    StepCollection<EigenproblemEnvironment> steps {};
    steps.add(DenseDiagonalization(number_of_requested_eigenpairs));

    return Algorithm<EigenproblemEnvironment>(steps);
}


}  // namespace EigenproblemSolver
}  // namespace GQCP
