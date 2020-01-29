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


#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"

#include <type_traits>


namespace GQCP {


/**
 *  An iteration step that calculates the subspace matrix.
 */
class SubspaceMatrixCalculation :
    public Step<EigenproblemEnvironment> {

public:

    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Calculate the subspace matrix, i.e. the projection of the matrix A onto the subspace spanned by the vectors in V.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment& environment) override {

        const auto& V = environment.V;   // the subspace of guess vectors
        const auto& VA = environment.VA;   // VA = A * V (implicitly calculated through the matrix-vector product)

        std::cout << "VA: " << std::endl << VA << std::endl << std::endl;

        environment.S = V.transpose() * VA;  // the "subspace matrix": the projection of the matrix A onto the subspace spanned by the vectors in V
        std::cout << "S: " << std::endl << environment.S << std::endl << std::endl;
    }
};


}  // namespace GQCP
