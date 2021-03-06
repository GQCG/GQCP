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
#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"
#include "QCModel/HF/UHF.hpp"


namespace GQCP {


/**
 *  An iteration step that calculates the current alpha and beta density matrices (expressed in the scalar/AO basis) from the current coefficient matrices.
 * 
 *  @tparam _Scalar              The scalar type used to represent the expansion coefficient/elements of the transformation matrix: real or complex.
 */
template <typename _Scalar>
class UHFDensityMatrixCalculation:
    public Step<UHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = UHFSCFEnvironment<Scalar>;


public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return A textual description of this algorithmic step.
     */
    std::string description() const override {
        return "Calculate the current UHF alpha and beta density matrices (in the AO basis) and place them in the environment.";
    }


    /**
     *  Calculate the current UHF alpha and beta density matrices (in the AO basis) and place them in the environment.
     * 
     *  @param environment              The environment that acts as a sort of calculation space.
     */
    void execute(Environment& environment) override {

        const auto& C = environment.coefficient_matrices.back();  // The most recent alpha and beta coefficient matrix.

        const auto D = QCModel::UHF<Scalar>::calculateScalarBasis1DM(C, environment.N.alpha(), environment.N.beta());

        environment.density_matrices.push_back(D);
    }
};


}  // namespace GQCP
