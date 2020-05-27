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
 *  An iteration step that calculates the alpha- and beta- error matrices from the Fock and density matrices (expressed in the scalar/AO basis).
 * 
 *  @tparam _Scalar              the scalar type used to represent the expansion coefficient/elements of the transformation matrix
 */
template <typename _Scalar>
class UHFErrorCalculation:
    public Step<UHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = UHFSCFEnvironment<Scalar>;


public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Calculate the current alpha- and beta- error vectors and add them to the environment.";
    }


    /**
     *  Calculate the current alpha- and beta- error vectors and add them to the environment.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        // Read F, D and S from the environment.
        const auto& S = environment.S;

        const auto& F_alpha = environment.fock_matrices_alpha.back();
        const auto& F_beta = environment.fock_matrices_beta.back();

        const auto& D_alpha = environment.density_matrices_alpha.back();
        const auto& D_beta = environment.density_matrices_beta.back();


        // Calculate the errors and write them to the environment (as a vector).
        const auto error_matrix_alpha = QCModel::UHF<Scalar>::calculateError(F_alpha, D_alpha, S);
        const auto error_matrix_beta = QCModel::UHF<Scalar>::calculateError(F_beta, D_beta, S);

        environment.error_vectors_alpha.push_back(error_matrix_alpha.pairWiseReduced());
        environment.error_vectors_beta.push_back(error_matrix_beta.pairWiseReduced());
    }
};


}  // namespace GQCP
