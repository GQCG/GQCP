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
 *  An iteration step that calculates the current density matrix (expressed in the scalar/AO basis) from the current coefficient matrix.
 * 
 *  @tparam _Scalar              the scalar type used to represent the expansion coefficient/elements of the transformation matrix
 */
template <typename _Scalar>
class RHFDensityMatrixCalculation :
    public IterationStep<RHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = RHFSCFEnvironment<Scalar>;


public:

    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Calculate the current RHF density matrix and place it in the environment
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {
        const auto& C = environment.coefficient_matrices.back();  // the most recent coefficient matrix
        const auto D = QCModel::RHF<double>::calculateScalarBasis1RDM(C, environment.N);
        environment.density_matrices.push_back(D);
    }
};


}  // namespace GQCP
