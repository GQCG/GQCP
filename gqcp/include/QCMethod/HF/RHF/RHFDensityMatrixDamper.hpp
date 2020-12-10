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
#include "Mathematical/Optimization/Accelerator/ConstantDamper.hpp"
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"
#include "QCModel/HF/RHF.hpp"


namespace GQCP {


/**
 *  An iteration step that accelerates the density matrix (expressed in the scalar/AO basis) based on a constant damping accelerator.
 * 
 *  @tparam _Scalar              The scalar type used to represent the expansion coefficient/elements of the transformation matrix: real or complex.
 */
template <typename _Scalar>
class RHFDensityMatrixDamper:
    public Step<RHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = RHFSCFEnvironment<Scalar>;


private:
    ConstantDamper damper;  // The damping accelerator.


public:
    /*
     *  CONSTRUCTORS 
     */

    /**
     *  @param alpha            The damping factor.
     */
    RHFDensityMatrixDamper(const double alpha) :
        damper {alpha} {}


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return A textual description of this algorithmic step.
     */
    std::string description() const override {
        return "Replace the most recent density matrix with an accelerated one.";
    }


    /**
     *  Replace the most recent density matrix with an accelerated one.
     * 
     *  @param environment              The environment that acts as a sort of calculation space.
     */
    void execute(Environment& environment) override {

        if (environment.density_matrices.size() < 2) {
            return;  // No acceleration is possible.
        }

        // Get the two most recent density matrices and produce an accelerated density matrix.
        const auto second_to_last_it = environment.density_matrices.end() - 2;
        const auto& D_previous = *second_to_last_it;  // Dereference the iterator.
        const auto D_current = environment.density_matrices.back();

        const auto D_accelerated = this->damper.accelerate(D_current, D_previous);
        environment.density_matrices.pop_back();  // We will replace the most recent density matrix with the accelerated one, so remove the most recent one.
        environment.density_matrices.push_back(D_accelerated);
    }
};


}  // namespace GQCP
