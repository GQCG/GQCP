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


#include "Mathematical/Algorithm/IterationStep.hpp"
#include "QCMethod/HF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/RHF.hpp"


namespace GQCP {


/**
 *  An iteration step that calculates the error matrix from the Fock and density matrices (expressed in the scalar/AO basis).
 * 
 *  @tparam _Scalar              the scalar type used to represent the expansion coefficient/elements of the transformation matrix
 */
template <typename _Scalar>
class RHFErrorCalculation :
    public IterationStep<RHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = RHFSCFEnvironment<Scalar>;


public:

    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Calculate the current error vector and add it to the environment.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        // Read F, D and S from the environment
        const auto& D = environment.density_matrices.back();
        const auto& S = environment.S;
        const auto& F = environment.fock_matrices.back();

        // Calculate the error and write it to the environment (as a vector)
        const auto error_matrix = calculateRHFError(F, D, S);
        environment.error_vectors.push_back(error_matrix.pairWiseReduce());
    }
};


}  // namespace GQCP
