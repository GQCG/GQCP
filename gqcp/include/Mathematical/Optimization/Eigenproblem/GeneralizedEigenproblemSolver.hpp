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

#pragma once


#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Optimization/Eigenproblem/GeneralizedDenseDiagonalization.hpp"


namespace GQCP {
namespace GeneralizedEigenproblemSolver {


/**
 *  @return An algorithm that can diagonalize a dense matrix.
 *
 *  @tparam Scalar          The scalar type of matrix elements: real or complex.
 */
template <typename Scalar>
Algorithm<GeneralizedEigenproblemEnvironment<Scalar>> Dense() {

    // Use Eigen's routines to diagonalize a matrix.
    StepCollection<GeneralizedEigenproblemEnvironment<Scalar>> steps {};
    steps.add(GeneralizedDenseDiagonalization<Scalar>());

    return Algorithm<GeneralizedEigenproblemEnvironment<Scalar>>(steps);
}


}  // namespace GeneralEigenproblemSolver
}  // namespace GQCP
