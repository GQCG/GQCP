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


#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"

#include <type_traits>


namespace GQCP {


/**
 *  An iteration step that calculates the subspace matrix, i.e. the projection of the matrix A onto the subspace spanned by the vectors in V.
 */
class SubspaceMatrixCalculation:
    public Step<EigenproblemEnvironment<double>> {

public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Calculate the subspace matrix, i.e. the projection of the matrix A onto the subspace spanned by the vectors in V.";
    }


    /**
     *  Calculate the subspace matrix, i.e. the projection of the matrix A onto the subspace spanned by the vectors in V.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment<double>& environment) override {

        const auto& V = environment.V;    // the subspace of guess vectors
        const auto& VA = environment.VA;  // VA = A * V (implicitly calculated through the matrix-vector product)

        environment.S = V.transpose() * VA;  // the "subspace matrix": the projection of the matrix A onto the subspace spanned by the vectors in V
    }
};


}  // namespace GQCP
