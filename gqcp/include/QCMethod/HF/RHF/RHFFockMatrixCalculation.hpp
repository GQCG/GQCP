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
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"
#include "QCModel/HF/RHF.hpp"


namespace GQCP {


/**
 *  An iteration step that calculates the current Fock matrix (expressed in the scalar/AO basis) from the current density matrix.
 * 
 *  @tparam _Scalar              The scalar type used to represent the expansion coefficient/elements of the transformation matrix: real or complex.
 */
template <typename _Scalar>
class RHFFockMatrixCalculation:
    public Step<RHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = RHFSCFEnvironment<Scalar>;


public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return A textual description of this algorithmic step.
     */
    std::string description() const override {
        return "Calculate the current RHF Fock matrix (expressed in the scalar/AO basis) and place it in the environment.";
    }


    /**
     *  Calculate the current RHF Fock matrix (expressed in the scalar/AO basis) and place it in the environment.
     * 
     *  @param environment              The environment that acts as a sort of calculation space.
     */
    void execute(Environment& environment) override {
        const auto& D = environment.density_matrices.back();  // The most recent density matrix.
        const auto F = QCModel::RHF<Scalar>::calculateScalarBasisFockMatrix(D, environment.sq_hamiltonian);
        environment.fock_matrices.push_back(F.parameters());
    }
};


}  // namespace GQCP
